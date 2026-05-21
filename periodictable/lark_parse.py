from typing import cast

import lark

from .core import PeriodicTable, Element, Atom, Isotope
from .core import default_table
from .formulas import (
    Formula, Structure,
    _mix_by_weight_pairs, _mix_by_volume_pairs,
    pretty as pretty_formula
)
from .util import from_subscript, from_superscript

# TODO: valence belongs to a group rather than element

# TODO: Parser can't handle meters as 'm' because it conflicts with the milli prefix
LENGTH_UNITS = {'nm': 1e-9, 'um': 1e-6, 'μm': 1e-6, 'mm': 1e-3, 'cm': 1e-2, 'Ang': 1e-10, 'Å': 1e-10}
MASS_UNITS = {'ng': 1e-9, 'ug': 1e-6, 'mg': 1e-3, 'g': 1e+0, 'kg': 1e+3}
VOLUME_UNITS = {'nL': 1e-9, 'uL': 1e-6, 'mL': 1e-3, 'L': 1e+0}

# TODO: use grammar string directly in the sphinx/guide/formula_grammar.rst
grammar = """
start      : SPACE? formula SPACE?  # strip blank space from start and end
formula    : compound | mixture

# Mixture definitions:  quantity compound // quantity compound // quantity compound
# Activation only cares about total mass, so you can freely mix masses and volumes if
# you have the density for each component. For scattering you need the density of the
# mixture. When this is different from the mixture of densities use (mixture)@density.
# For thin film samples, allow stacking of layers with the thickness of each layer.
# With density for each layer the relative quantities of each element in the stack can
# be calculated. Convert to mass by multiplying by thickness (cm) and area (cm²).

mixture    : byamount | byvolume | byweight | layers
byamount   : quantity compound (MIX quantity compound)*
byvolume   : volumepct compound (MIX percentage compound)* MIX compound
byweight   : weightpct compound (MIX percentage compound)* MIX compound
layers     : thickness compound (MIX thickness compound)*
quantity   : NUMBER SPACE? (MASS | VOLUME) SPACE
weightpct  : NUMBER SPACE? WEIGHTPCT SPACE
volumepct  : NUMBER SPACE? VOLUMEPCT SPACE
thickness  : NUMBER SPACE? LENGTH SPACE
percentage : NUMBER SPACE? "%" SPACE  # Allows "3 % "

# Composite: number group ... @density where group is El count El count ...
# Density applies to the entire composite, such as "NaCl + 29.2H2O @ 1.07n"
# For the density of a mixture you need parentheses: "(10 wt% NaCl // H2O)@1.07n"
# FASTA sequences: (rna|dna|aa) : SEQUENCE @ density
# Note: optional `[token]` leaves a None placeholder in the tree, unlike `token?`
compound   : (composite | fasta) [density]
fasta      : FASTA ":" SEQUENCE
composite  : [NUMBER] group (SEPARATOR [NUMBER] group)*
group      : ((atom | isoatom | "(" formula ")") [COUNT])+
atom       : SYMBOL [isotope] [valence]
isoatom    : SUPERINT SYMBOL [valence]    # For example ²H for deuterium
isotope    : "[" INTEGER "]"
valence    : "{" [INTEGER] CHARGE "}" | [SUPERINT] SUPERCHARGE
density    : SPACE? "@" SPACE? DENSITY [DENSITYMODE]

# Tokens
#FASTA     : /dna|rna|aa/  # Sequence type is limited to these values but ...
FASTA      : /[a-z]+/      # "str:sequence" syntax allows better error reporting
SEQUENCE   : /[-A-Z *]+/
# could list all elements, but better error reporting if element symbol lookup fails
SYMBOL     : /[A-Z][a-z]*/
CHARGE     : /[+]+|[-]+/  # allow valence using {++} or {--}
DENSITY    : NUMBER  # using alias DENSITY for number for better error reporting
DENSITYMODE: /[ni]/       # n=natural density, i=isotopic density
MIX        : SPACE? "//" SPACE?
WEIGHTPCT  : /%w((eigh)?t)?/ | /w((eigh)?t)?%/ | /%m(ass)?/ | /m(ass)?%/
VOLUMEPCT  : /%v(ol(ume)?)?/ | /v(ol(ume)?)?%/
MASS       : "kg" | "g" | "mg" | "ug" | "μg" | "ng"
VOLUME     : "L" | "mL" | "uL" | "μL" | "nL"
LENGTH     : "cm" | "mm" | "um" | "μm" | "nm" | "Ang" | "Å"
COUNT      : NUMBER | SUBNUM  # atom counts can be normal numbers or unicode subscripts

SEPARATOR  : SPACE? /[+•·]/ SPACE? | SPACE
SPACE      : /[ \\t\\n\\r]+/
NUMBER     : INTEGER | FRACTION
INTEGER    : /[1-9][0-9]*/
FRACTION   : /([1-9][0-9]*|0)?[.][0-9]*/  # allow all floats?
SUBNUM     : SUBINT | SUBFRAC
SUBINT     : /(\u2080|[\u2081-\u2089][\u2080-\u2089]*)/
SUBFRAC    : /(\u2080|[\u2081-\u2089][\u2080-\u2089]*|)([.][\u2080-\u2089]*)/
SUPERINT   : /(\u2070|[\u00B9\u00B2\u00B3\u2074-\u2079][\u2070\u00B9\u00B2\u00B3\u2074-\u2079]*)/
SUPERCHARGE: /\u207A+|\u207B+/  # Allow Ca++ and Cl- using superscript + and -
"""

# propagate_positions saves start_pos and end_pos for each rule as well as each terminal.
formula_parser = lark.Lark(grammar, propagate_positions=True)

def int_or_float(s):
    f = float(s)
    i = int(f)
    return i if i == f else f

class StripPunctuation(lark.Transformer):
    """
    Token stripper visitor class.

    This is done separately from the formula composer so that we can show the cleaned tree
    before debugging the conversion.

    Unnamed punctuation characters []{}():% and units (kg, mL, nm, ...) which are represented
    as quoted strings in the grammar have no associated token.

    Note: could get the same effect by renaming the unused terminals with leading underscore,
    but that makes the grammar harder to read.
    """
    def SEPARATOR(self, _):
        """Strip token for molecular fragment separator (+ or center dot or spaces)."""
        return lark.Discard
    def MIX(self, _):
        """Strip token for mixture separator //."""
        return lark.Discard
    def SPACE(self, _):
        """Strip token for (usually optional) spaces."""
        return lark.Discard
    def WEIGHTPCT(self, _):
        """Strip token for wt% mixture indicator."""
        return lark.Discard
    def VOLUMEPCT(self, _):
        """Strip token for vol% mixture indicator."""
        return lark.Discard

class ConvertTokens(lark.Transformer):
    """
    Syntax tree to formula conversion class.
    """
    def __init__(self, text, table=None):
        """
        *text* is the original formula string.

        *table* is an optional alternative periodic table.
        """
        self._context = text
        self._table = default_table(table)

    def VOLUME(self, token: lark.Token) -> tuple[str, str]:
        """
        Convert VOLUME terminal ('volume', unit) pair.

        Unit is a volume unit, such as mL or uL for microlitres.
        """
        return 'volume', token.value
    def MASS(self, token: lark.Token) -> tuple[str, str]:
        """
        Convert MASS terminal to ('mass', unit) pair.

        Unit is a mass unit, such as g or mg.
        """
        return 'mass', token.value
    def LENGTH(self, token: lark.Token) -> tuple[str, str]:
        """
        Convert LENGTH terminal to ('length', unit) pair.

        Unit is a length unit, such as cm or nm.
        """
        return 'length', token.value
    def NUMBER(self, token: lark.Token) -> int|float:
        """
        Convert string to float or integer.

        Numbers are used for quantities and percentages in mixtures, and for multiplier
        counts to molecule fragments.
        """
        return int_or_float(token.value)
    DENSITY = NUMBER  # We've aliased DENSITY and NUMBER in the grammar
    def INTEGER(self, token: lark.Token) -> int:
        """
        Convert string to float or integer
        """
        return int(token.value)
    def COUNT(self, token: lark.Token) -> int|float:
        """
        Return the count value for a group component.

        Count is specified after the symbol, either as an ASCII number or using subscript digits.
        The period separator for fractional counts uses ASCII in both cases (there is no subscript
        period charcter available). If the count is fractional return it as a float, otherwise
        return it as an integer.
        """
        return int_or_float(from_subscript(token.value))
    def SUPERINT(self, token) -> int:
        """
        Return the integer value of a sequence of superscript digits.

        This is used to specify the valence or to specify the isotope.
        """
        return int(from_superscript(token.value))
    def DENSITYMODE(self, token) -> str:
        """
        Return the value of the DENSITYMODE token, either "n" or "i". If no mode is specified
        then a token value of None will be given to the density rule.
        """
        return token.value
    def CHARGE(self, token) -> int:
        """
        Return a sequence of plus and minus characters. By grammar rules they must all have
        the same sign.

        This is used in the valence rule to specify the charge for the atom.
        """
        return token.value
    def SUPERCHARGE(self, token) -> str:
        """
        Convert sequence of superscript plus and minus characters to ASCII plus and minus.

        This is used in the valence rule to specify the charge for the atom.
        """
        return from_superscript(token.value)
    def SYMBOL(self, token) -> Element:
        """
        Look up the element in the periodic table and return it.

        Raise ValueError if the element doesn't exist.
        """
        try:
            return self._table.symbol(token.value)
        except Exception:
            raise ValueError(f"Element {token.value} doesn't exist")
    def FASTA(self, token) -> str:
        """
        Return the token value as the fasta sequence type: "dna", "rna" or "aa".
        """
        return token.value
    def SEQUENCE(self, token) -> str:
        """
        Return the token value as the fasta sequence string.
        """
        return token.value
    def isotope(self, tokens) -> int:
        """
        Return the isotope number for the atom.

        Transform: [isotope] => isotope
        """
        return tokens[0]
    def valence(self, tokens) -> int:
        """
        Return valence from number and sign.

        Valence is either a number followed by plus or minus, or a sequence of plus
        or minus. If the number was specified it will already have been converted
        to a value, otherwise use the length of the charge string as the value.

        The valence can be given using superscript or regular ASCII number and sign
        symbols. If ASCII then they need to be wrapped in braces such as Ca{2+}. The
        token transform handles the conversion from superscript to ASCII characters
        and the conversion from string to number.

        Raise ValueError if a number was supplied along with multiple charge symbols.

        Transform: [number|None, 'charge'] => valence

        Example: ['{1+}'] => [1, '+'] = Ca.ion[1]

        Example: Ca{++} => [None, '++'] = Ca.ion[2]

        Example: Ca{3--} => ValueError
        """
        # print("in valence with", tokens)
        value, charge = tokens
        if value is None:
            value = len(charge)
        elif value and len(charge) > 1:
            raise ValueError(f"Use {value}{charge[0]} instead of {value}{charge} for valence")
        valence = value if charge[0] == '+' else -value
        return valence
    def atom(self, tokens) -> Atom:
        """
        Returns an atom from the periodic table.

        Usually this will use elements from the default table, but if an alternate table is
        provided to the ConvertTokens constructor then that will be used to retrieve the element
        from the symbol.

        Raises an error if the symbol does not exist, does not have that isotope or doesn't
        allow that valence.

        Transform: ['symbol', isotope|None, valence|None] => atom

        Example: ['H', 1, 1] => H[1]{+}

        Example: ['Ca', None, 2] => Ca{2+}
        """
        #print("atom", tokens)
        el, iso, ion = tokens
        if iso and ion:
            atom = el[iso].ion[ion]
        elif iso:
            atom = el[iso]
        elif ion:
            atom = el.ion[ion]
        else:
            atom = el
        #print(f"atom {tokens} => {atom}")
        return atom

    def isoatom(self, tokens) -> Atom:
        """
        Returns an isotope from the periodic table.

        Usually this will use elements from the default table, but if an alternate table is
        provided to the ConvertTokens constructor then that will be used to retrieve the element
        from the symbol.

        Raises an error if the symbol does not exist, does not have that isotope or doesn't
        allow that valence.

        Transform: [isotope, 'symbol', valence|None] => atom

        Example ²H⁺: [2, 'H', 1] => D{+}
        """
        # print("isoatom", tokens)
        iso, el, ion = tokens
        atom = el[iso].ion[ion] if ion else el[iso]
        # print(f"isoatom {tokens} => {atom}")
        return atom


    def group(self, tokens) -> Structure:
        """
        Returns a sequence of (count, item) pairs, where item is an atom or a nested formula.
        Missing counts default to 1.

        Transform: [atom|formula, count|None, ...] => ((count, atom|formula), ...)

        Example CaCO3: [Ca, None, C, None, O, 3]
        => ((1, Ca), (1, C), (3, O))
        """
        # print("group tokens", tokens)
        tokens = [1 if value is None else value for value in tokens]
        pairs = tuple((count, item) for item, count in zip(tokens[::2], tokens[1::2]))
        # print("group output", pairs)
        return pairs

    def composite(self, tokens) -> Structure:
        """
        Returns a sequence of (number, group) pairs. Each group is a sequence of (count, item)
        pairs, where item is an atom or a nested formula. Missing numbers default to 1.

        Transform: [number|None, group, ...] => ((number, group), ...) | ((count, atom), ...)

        Example CaCO3 6H2O: [None, ((1, Ca), (1, C), (3, O)), 6, ((2, H), (1, O))]
        => ((1, ((1, Ca), (1, C), (3, O))), (6, ((2, H), (1, O))))

        Example CaCO3(H2O)6: [None, ((1, Ca), (1, C), (3, O), (6, formula('H2O')))]
        => ((1, ((1, Ca), (1, C), (3, O), (6, formula('H2O')))),)

        Example CaCO3 (H2O)6: [None, ((1, Ca), (1, C), (3, O)), None, ((6, formula('H2O')),)]
        => ((1, ((1, Ca), (1, C), (3, O))), (1, ((6, formula('H2O')),)))
        """
        # print("composite tokens", tokens)
        numbers = [1 if v is None else v for v in tokens[::2]]
        groups = tokens[1::2]
        pairs = tuple((number, group) for number, group in zip(numbers, groups))
        # print("composite output", pairs)
        return pairs

    def fasta(self, tokens) -> Structure:
        r"""
        Returns the formula corresponding to the FASTA sequence, with the natural
        density set. Labile hydrogen use H[1] in the formula.

        The extra level of nesting in the return value is so that the fasta structure
        is like a composite with a single group containing a nested formula.

        Transform: [ 'aa|dna|rna', '[-A-Z \*]+' ] => (1, ((1, formula),))

        Example: dna:CAGT: ['dna', 'CAGT'] x=> ((1, ((1, formula('C39H37H[1]10N15O25P4')),)),)
        """
        # TODO: fasta is ignoring table when parsing
        # TODO: avoid circular imports
        # TODO: support other biochemicals (carbohydrate residues, lipids)
        from periodictable.fasta import CODE_TABLES, Sequence

        # print("fasta input", tokens)
        seq_type, seq = tokens
        if seq_type not in CODE_TABLES:
            raise ValueError(f"Invalid fasta sequence type '{seq_type}:'")
        seq = Sequence(name="seq", sequence=seq, type=seq_type)
        pairs = ((1, seq.labile_formula),)
        composite = ((1, pairs),)
        # print("fasta output", composite)
        # return tuple[tuple[int, tuple[tuple[int, Formula]]]] as Structure
        return cast(Structure, composite)

    def density(self, tokens) -> tuple[str, float, str]:
        """
        Returns a density tuple from the @density construct. Density mode 'n' for
        natural or 'i' for isotopic defaults to isotopic. That is, D2O@1.11 is the
        isotopic density of D2O, not the natural density of H2O with conversion to
        the heavier deutrium isotope.

        Transform: [value, mode|None] => ('density', value, mode)

        Example @1.11: [1.11, None] => ('density', 1.11, 'i')

        Example @1.11i: [1.11, 'i'] => ('density', 1.11, 'i')

        Example @1n: [1, 'n'] => ('density', 1, 'n')
        """
        value = tokens[0]
        mode = 'i' if not tokens[1] else tokens[1]
        return 'density', value, mode

    def compound(self, tokens) -> Formula:
        """
        Returns the formula for the compound, with optional density set.

        Density is ('density', value, mode) or None, where mode is 'i' for isotopic density
        or 'n' for natural density.

        The compound may come from a FASTA spec, such as dna:CAGT or from a composite, such
        as CaCO3+6H2O. The composite may include an embedded formula, such as CaCO3(H2O)6.
        In any case, the resulting material token will be a sequence of (multiplier, group)
        pairs, where each group is a sequence of (count, item) pairs. Each item may be an
        atom or a formula. The fasta transform returns a single group with a single item.
        As a nested sequence this is ((1, ((1, formula), ...)), ...), with nothing in the
        ellipses.

        Transform: [((number, group), ...), ('density', value, mode)|None] => formula

        Example NaCl@2.16i: [((1, ((1, Na), (1, Cl))),), ('density', 2.16, 'i')] => NaCl@2.16i

        Example dna:CAGT: [((1, ((1, formula('C39H37H[1]10N15O25P4')),)),), None] => C39H37H[1]10N15O25P4@1.69n

        Example CaCO3 6H2O: [((1, ((1, Ca), (1, C), (3, O))), (6, ((2, H), (1, O)))), None] => CaCO3(H2O)6

        Example CaCO3(H2O)6: [((1, ((1, Ca), (1, C), (3, O), (6, formula('H2O')))),), None] => CaCO3(H2O)6
        """
        # print("compound tokens", tokens)
        components, density_tuple = tokens
        if density_tuple is None:
            density, density_mode = None, 'i'
        else:
            _, density, density_mode = density_tuple

        # If a singleton formula with no density override then return it
        # That is, [(1, ((1, formula),)), None] => formula
        if density is None and len(components) == 1:
            number, group = components[0]
            if len(group) == 1 and number == 1:
                count, item = group[0]
                if count == 1 and isinstance(item, Formula):
                    # print("isolated formula with no density override")
                    return item

        # Not an isolated formula, so expand formulas within the groups.
        # That is, [..., (number, (..., (count, formula), ...)), ...]
        # becomes [..., (number, (..., (count, formula.structure), ...)), ...]
        def expand_formula(group):
            return tuple((count, getattr(item, 'structure', item)) for count, item in group)
        components = tuple((number, expand_formula(group)) for number, group in components)

        # If it is a singleton group then use its structure as the formula structure.
        if len(components) == 1 and components[0][0] == 1:
            structure = components[0][1]
        else:
            structure = components

        # Build the formula and assign density if available.
        # print("compound structure", structure)
        formula = Formula(structure=structure)
        if density is not None:
            if density_mode == 'n':
                formula.natural_density = density
            else:
                formula.density = density

        # print(f"compound output {formula} @ {formula.density}")
        return formula

    def weightpct(self, tokens) -> float:
        """
        Returns the percentage. The value has already be converted to a number.

        Used as the first percentage of a mix by weight mixture.

        Transform: [percent] => percent

        Example for "3 wt%": [3] => 3
        """
        return tokens[0]

    def volumepct(self, tokens) -> float:
        """
        Returns the percentage. The value has already be converted to a number.

        Used as the first percentage of a mix by volume mixture.

        Transform: [percent] => percent

        Example for "3 vol%": [3] => 3
        """
        return tokens[0]

    def percentage(self, tokens) -> float:
        """
        Returns the percentage. The value has already be converted to a number.

        Transform: [percent] => percent

        Example for " 3 % ": [3] => 3
        """
        return tokens[0]

    def byweight(self, tokens) -> Formula:
        """
        Returns mixture by wt% of the various components in the system.

        Raises ValueError if total exceeds 100%.

        Transform: [percent, formula, ..., percent, formula, formula] => formula

        Example: [76.95, D2O, H2O] => (D2O)3H2O
        """
        # TODO: structure not preserved in mixtures
        total = sum(tokens[:-1:2])
        if total > 100:
            raise ValueError(f"Total weight {total}% is more than 100% in wt% mixture")
        pairs = [(compound, percent) for percent, compound in zip(tokens[:-1:2], tokens[1:-1:2])]
        pairs.append((tokens[-1], 100-total))
        # return 'byweight', [*pairs, last_pair]
        formula = _mix_by_weight_pairs(pairs)
        # print(f"byweight => {formula} @ {formula.density}")
        return formula

    def byvolume(self, tokens) -> Formula:
        """
        Returns mixture by vol% of the various components in the system. Volumes are converted
        to mass using density.

        Raises ValueError if the density is missing from a component formula.
        Raises ValueError if total exceeds 100%.

        Transform: [percent, formula, ..., percent, formula, formula] => formula

        Example: [75.0, D2O@1n, H2O@1n] => (D2O)3H2O
        """
        # print("by volume", tokens)
        total = sum(tokens[:-1:2])
        if total > 100:
            raise ValueError(f"Total volume {total}% is more than 100% in vol% mixture")
        pairs = [(compound, percent) for percent, compound in zip(tokens[:-1:2], tokens[1:-1:2])]
        pairs.append((tokens[-1], 100-total))
        # print("byvolume pairs", pairs)
        # print("byvolume density", [f.density for f, p in pairs])
        #return 'byvolume', pairs
        formula = _mix_by_volume_pairs(pairs)
        return formula

    def byamount(self, tokens) -> Formula:
        """
        Returns mixture by mass of the various components in the system. Volumes are converted
        to mass using density.

        Raises ValueError if the density is missing from a component formula.

        Transform: [quantity, formula, ...] => formula

        Example: [('mass', 5.07, 'g'), NaCl@2.16, ('volume', 50, 'mL'), H2O@1n] => NaCl(H2O)32
        """
        # print("byamount", tokens)
        def find_value(quantity, formula):
            qtype, value, units = quantity
            if qtype == 'volume':
                if formula.density is None:
                    raise ValueError(f"Need the mass density of {formula}")
                mass = value * VOLUME_UNITS[units] * 1000.0 * formula.density
            else:
                mass = value * MASS_UNITS[units]
            return mass
        values = [find_value(q, f) for q, f in zip(tokens[::2], tokens[1::2])]
        total = sum(values)
        percent = [(m/total)*100 for m in values]
        formula = _mix_by_weight_pairs(zip(tokens[1::2], percent))
        formula.total_mass = total
        return formula

    def layers(self, tokens) -> Formula:
        """
        Returns the mixture by volume of the various layers in the system.

        Raises ValueError if the density is missing from a component formula.

        Sets formula.thickness to the sum of the layer thicknesses.

        Transform: [quantity, formula, ...] => formula

        Example: [('length', 10.006, 'nm'), Ni, ('length', 3, 'mm'), Si] => NiSi164000
        """
        # # Sanity check: make sure all units are length units. This won't happen
        # # because the parser only accepts proper formulas.
        # assert all(units in LENGTH_UNITS for dim, value, units in tokens[::2])
        values = [value*LENGTH_UNITS[units] for dim, value, units in tokens[::2]]
        total = sum(values)
        percent = [(m/total)*100 for m in values]
        formula = _mix_by_volume_pairs(zip(tokens[1::2], percent))
        formula.thickness = 100*total # convert meters to centimeters for cgs units
        return formula

    def mixture(self, tokens) -> Formula:
        """
        Returns the formula representing the mixture, either byweight, byvolume, byamount or layers

        Transform: [formula] => formula
        """
        return tokens[0]

    def formula(self, tokens) -> Formula:
        """
        Return the formula representing the compound or mixture.

        Transform:  [formula] => formula
        """
        return tokens[0]

    def thickness(self, tokens) -> tuple[str, float, str]:
        """
        Returns (dimension, value, unit) with dimension equal 'length'

        Transform: [value, ('length', unit)] => ('length', value, unit)

        Example: [5, ('length', 'nm')] => ('length', 5, 'nm')
        """
        value, (dim, units) = tokens
        return dim, value, units

    def quantity(self, tokens) -> tuple[str, float, str]:
        """
        Returns (dimension, value, unit) with dimension equal 'mass' or 'volume'

        Transform: [value, (dimension, unit)] => (dimension, value, unit)

        Example: [5, ('mass', 'g')] => ('mass', 5, 'g')
        """
        value, (dim, units) = tokens
        return dim, value, units

    def start(self, tokens) -> Formula:
        """
        Return the final formula, with the original text attached.

        Sets formula.name to the parser input string before returning.

        Transform: [formula] => formula
        """
        formula = tokens[0]
        # Remember the string which was parsed
        formula.name = self._context
        return formula

# TODO: if the next character is ":" then report error as bad fasta sequence type
def _allowed(allowed):
    # * SPACE, SEPARATOR: Generally ignored
    # * LPAR occurs whereever a symbol could be expected, so skip it
    # * COLON: If asking then it probably thinks it is looking for a fasta sequence, but
    # instead it should be looking for an element, so replace COLON with SYMBOL.
    # * AT: Looking for @DENSITY
    # * LPAR, RPAR: "(" and ")" are more readable
    # * LSQB: end of element, looking for isotope, so skip
    # * LBRACE, SUPERINT, SUPERCHARGE: end of element, looking for valence, so skip
    skip = set("SPACE SEPARATOR LPAR LSQB LBRACE SUPERINT SUPERCHARGE".split())
    # TODO: use order of elements in subst to sort the allowed list (currently alphabetical)
    subst = dict(
        NUMBER="NUMBER", # start of compound or start of mixture
        #FASTA="[dna|rna|aa]:SEQ",
        FASTA="aa:SEQ",
        COLON=":SEQ",
        #COLON="aa:SEQ",
        SEQUENCE="aa:SEQ",
        SEPARATOR="+", # generic group separator in composite
        SPACE="SPACE",
        SYMBOL="SYMBOL",
        CHARGE="CHARGE[+-]",
        LPAR='(',
        RPAR=')',
        LSQB='[',
        RSQB=']',
        LBRACE='{', # equivalent to SUPERINT and SUPERCHARGE
        RBRACE='}',
        VOLUMEPCT="vol%",
        WEIGHTPCT="wt%",
        MASS="UNIT[mg]",
        VOLUME="UNIT[mL]",
        LENGTH="UNIT[mm]",
        PERCENT="%",
        # I don't think all three of these can be concurrently allowed so no need to
        # deduplicate. Moot since the set operation happens again after substition below.
        AT="@DENSITY[ni]", # only the @ is expected, but better for doc
        DENSITY="@DENSITY[ni]", # only the number is expected, but better for doc
        DENSITYMODE="@DENSITY[ni]", # only the [ni] is expected, but better for doc
        MIX="//",
        # SUBNUM SUBINT SUBFRAC covered by COUNT
        # INTEGER and FRACTION covered by NUMBER
        # SUPERINT SUPERCHARGE LSQB LBRACE coexist with COUNT so stripped
        SUPERCHARGE="SUPERSCRIPT[+-]", # If you see a superscript number then you need a sign
        )
    stripped = set(s for s in allowed if s not in skip)
    if not stripped:
        stripped = allowed
    # Perform substitution for document strings
    stripped = set(subst.get(s, s) for s in stripped)
    if len(stripped) > 1:
        message = f"one of {' '.join(sorted(stripped))}"
    elif stripped:
        message = [*stripped][0]
    else:
        # This occurs when the middle part of percent mixtures have no percentage.
        # We could look for '//' in the string to report a better error message.
        message = "end of formula"
    return message

def parse_formula(formula_str: str, table: PeriodicTable|None=None) -> Formula:
    """
    Parse a chemical formula, returning a structure with elements from the
    given periodic table.
    """
    cleanup = StripPunctuation()
    convert = ConvertTokens(formula_str, table=table)
    try:
        tree = formula_parser.parse(formula_str)
    except lark.exceptions.UnexpectedCharacters as exc:
        #import pprint; pprint.pprint(exc.__dict__)
        context = exc.get_context(formula_str).rstrip()
        #context = exc._context.rstrip()
        message = f"Expected {_allowed(exc.allowed)} in\n{context}"
        raise ValueError(message)
    except lark.exceptions.UnexpectedEOF as exc:
        # import pprint; pprint.pprint(exc.__dict__)
        context = exc.get_context(formula_str).rstrip()
        message = f"Expected {_allowed(exc.expected)} in\n{context}"
        raise ValueError(message)
    except Exception as exc:
        # TODO: are other exceptions possible from the Earley parser?
        raise exc from None
    tree = cleanup.transform(tree)
    try:
        tree = convert.transform(tree)
    except lark.exceptions.VisitError as exc:
        # Unwind the VistorError exception capture and reraise the original exception
        # This requires that error messages in the transformer give enough context to
        # correct the error.
        raise exc.orig_exc from None
    return tree

# Error conditions are marked with '!' so the exception is ignored
# Lines marked ## fail on the existing parser
examples = """
! DNA:CAGT  # incorrect case for FASTA type not properly identified
! dna CAGT  # missing colon in FASTA
! O²  # SUPERCHARGE should be the only valid token here
! ₃H2O  # badly placed subscript
! // 3g Ca  # // is not a comment
! 3g Ca@ // 5g Si # missing density value
! Ca@i  # missing density value  ##
! Ca ⁺⁺  # extra space before valence
! Ca++  # missing braces in valence: the + is acting as SEPARATOR
! Ca2+  # missing braces in valence: the 2 is acting as COUNT and the + as SEPARATOR
! Ca{2}  # missing charge in valence
! 37 vol% H2O@1 / 5% D2O@1  # missing /
! 37 vol% H2O@1 /// 5% D2O@1  # extra /
! H2O@1h  # bad density mode
! 37 vol% NaCl@2.16 // H2O@1 // D2O@1  # percent missing in middle part
! 37 vol% H2O@1 // 5% D2O@1  # percent not allowed in last part
! 37 vol% H2O@1 // 5 vol% D2O@1  # only % in subsequent parts
! 37% H2O@1 // D2O@1  # missing vol% or wt%
! 37 val% H2O@1 // D2O@1  # bad spelling of vol%
! Fe[56O2 # bad isotope syntax
! Co[181]  # bad isotope
! Ca{2+O2  # bad valence syntax
! Co{17-}  # bad valence
! 3..5 mg NaCl
! 3.5 fm Si # bad units at the start; could be wt%/vol% or LENGTH, VOLUME, MASS 
! 3.5 mm Si // 2.5 nm SiO2 //
! 3.5 mm Si // 2.5 nm SiO2 // 35 mm cG
! ((Co) # mismatched LPAR
! Co)  # mismatched RPAR
! bad:CAGT  # bad sequence type
Co
dna:CAGT
(Co@5) ##
(((Co@5)@6)) ##
CaCO3
CaCO₃
CaCO3+6H2O
CaCO3 6H2O
CaCO3(H2O)6
CaCO3 (H2O)6
(Ca(CO3)((H2O)6))
CaCO₃·6H₂O  ##
DHO
!Ca{2++}  # bad valence string
Ca⁺⁺  # also Ca{2+}  ##
O²⁻   ##
H[1]
²H⁺    # D{+} ##
O²H⁻   # OD{-} ##
O²⁻H⁺  # O{2-}H{+} ##
O²⁻²H⁺ # O{2-}D{+} ##
H2O@1
D2O@1n
D2O @ 1.11  ##
D2O@1.11i
HO{1-}
H[1]{1-}O
H2SO4
C3H4H[1]NO@1.29n
78.2H2O[16] + 21.8H2O[18] @1n  # density applies to composite
dna:CAGT @1n  # fasta density override
50 wt% Co // Ti
33 wt% Co // 33% Fe // Ti
! 93 wt% Co // 33% Fe // Ti  # More than 100 wt%
! 93 vol% Co // 33% Fe // Ti  # More than 100 vol%
20 vol% (10 wt% NaCl@2.16 // H2O@1) // D2O@1n
NaCl(H2O)29.1966(D2O)122.794@1.10i
5g NaCl // 50mL H2O@1
5g NaCl@2.16 // 50mL H2O@1
! 5g NaCl // 50mL H2O   # Need density for H2O to convert volume to mass
(10 wt% NaCl // H2O)@1.07n # set density of a mixture
50 mL (45 mL H2O@1 // 5 g NaCl)@1.0707 // 20 mL D2O@1n
1 cm Si // 5 nm Cr // 10 nm Au
aa:RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV

! Bl2Oh   # Bad symbol
! 5 Mg NaCl // 50mL H2O@1  # Bad units
! 4 nm NaCl@2.17// 50 g Si  # Can't use mass in layer mixture

"""

def check():
    from periodictable.formulas import old_parser

    for line in examples.split('\n'):
        formula = line.split('#')[0]
        bad = line.startswith('!')
        if bad:
            formula = formula[1:]
        if formula:
            if bad:
                print(f"!!! {line[1:]}")
            else:
                print(f"*** {line}")
            try:
                # Toggle the following to test pyparsing vs lark
                tree = parse_formula(formula)
                #tree = old_parser(formula) if "##" not in line else "!!! pyparsing fails"
                density = getattr(tree, 'density', None)
                density_str = f" @ {density:.2f}" if density else ""
                mode = 'unicode' # unicode latex html plain
                # mode = 'plain'
                print(f" => {pretty_formula(tree, mode)}{density_str}")
                # print(f"    {getattr(tree, 'structure', None)}")
            except Exception as exc:
                if bad:
                    print(f"{exc}")
                else:
                    raise exc from None
            else:
                if '##' in line:
                    continue  # pyparsing should fail but doesn't
                if bad:
                    raise RuntimeError(f"Exception not raised for <{formula}>")

def main():
    import sys

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            formula = parse_formula(arg)
            mass = f" {formula.total_mass:.4g} g" if formula.total_mass else ""
            density = f"@{formula.density:.4g}" if formula.density else ""
            thickness = f" {10*formula.thickness:.4g} mm" if formula.thickness else ""
            print(f"{formula}{density}{mass}{thickness}")
    else:
        check()

if __name__ == "__main__":
    main()
