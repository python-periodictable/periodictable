import lark
import periodictable as pt
from periodictable.core import PeriodicTable
from periodictable.core import default_table
from periodictable.formulas import from_subscript, Formula, _mix_by_weight_pairs, _mix_by_volume_pairs
from periodictable.formulas import VOLUME_UNITS, MASS_UNITS, LENGTH_UNITS

grammar = """
start      : SPACE? formula SPACE? # strip blank space from start and end
formula    : compound | mixture

# Mixture definitions:  quantity compound // quantity compound // quantity compound
# Activation only cares about total mass, so you can freely mix masses and volumes if
# you have the density for each component. Scattering cares about density of the mixture,
# which in general is different from the mixture of densities.
# To convert layers to masses for activation estimates we need density. Also need to scale by
# area to convert density and thickness to mass. Assume unit area is cm^2, so for
# example "4 (5 nm Ni // 2 mm Si)" is a 4 cm^2 wafer of nickel on silicon. If you
# were to add a polymer you would need its density: "4 (20 nm C5H10@1.2

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

# Compound definition: number group ... @ density where group is El count El count ...
# FASTA sequences: (rna|dna|aa) : SEQUENCE @ density
# Density applies to the entire formula, such as "NaCl + 29.2H2O @ 1.07n"
# If you do this as a mixture you need parentheses: "(10 wt% NaCl // H2O)@1.07n"
# Note: `[token]` leaves a None placeholder in the tree, unlike `token?`
compound   : (composite | fasta) [density]
fasta      : FASTA ":" SEQUENCE
FASTA      : /dna|rna|aa/
SEQUENCE   : /[A-Z -*]+/
composite  : [NUMBER] group (SEPARATOR [NUMBER] group)*
group      : ((atom | "(" formula ")") [COUNT])+
atom       : SYMBOL [isotope] [charge]
# could list all elements, but better error reporting if element symbol lookup fails
SYMBOL     : /[A-Z][a-z]*/
isotope    : "[" INTEGER "]"
charge     : "{" [INTEGER] CHARGE "}" | [SUPERINT] SUPERCHARGE
density    : SPACE? "@" SPACE? NUMBER [DENSITYMODE]

# Tokens
CHARGE     : /[+]+|[-]+/  # allow charge using {++} or {--}
SUPERCHARGE: /\u207A+|\u207B+/  # Allow Ca++ and Cl- using superscript + and -
DENSITYMODE: /[ni]/
MIX        : SPACE? "//" SPACE?
# maybe drop "wt%" and "vol%"
WEIGHTPCT  : /%w((eigh)?t)?/ | /w((eigh)?t)?%/ | /%m(ass)?/ | /m(ass)?%/
VOLUMEPCT  : /%v(ol(ume)?)?/ | /v(ol(ume)?)?%/
MASS       : "kg" | "g" | "mg" | "ug" | "μg" | "ng"
VOLUME     : "L" | "mL" | "uL" | "μL" | "nL"
LENGTH     : "cm" | "mm" | "um" | "μm" | "nm"

SEPARATOR  : SPACE? /[+•·]/ SPACE? | SPACE
SPACE      : /[ \\t\\n\\r]+/
COUNT      : NUMBER | SUBNUM  # atom counts can be normal numbers or unicode subscripts
NUMBER     : INTEGER | FRACTION
INTEGER    : /[1-9][0-9]*/
FRACTION   : /([1-9][0-9]*|0)?[.][0-9]*/  # allow all floats?
SUBNUM     : SUBINT | SUBFRAC
SUBINT     : /(\u2080|[\u2081-\u2089][\u2080-\u2089]*)/
SUBFRAC    : /(\u2080|[\u2081-\u2089][\u2080-\u2089]*|)([.][\u2080-\u2089]*)/
SUPERINT   : /(\u2070|[\u00B9\u00B2\u00B3\u2074-\u2079][\u2070\u00B9\u00B2\u00B3\u2074-\u2079]*)/
"""

parser = lark.Lark(grammar)

def from_superscript(value: str) -> str:
    """
    Convert unicode superscript characters to normal characters. This allows us to parse,
    for example, Ca²⁺ as Ca{2+}.
    """
    codepoints = {
        '\u2070': '0', '\u00B9': '1', '\u00B2': '2', '\u00B3': '3',
        '\u2074': '4', '\u2075': '5', '\u2076': '6', '\u2077': '7',
        '\u2078': '8', '\u2079': '9', '\u207a': '+', '\u207b': '-',
        '\u207c': '=', '\u207d': '(', '\u207e': ')',

        '\u2071': 'i', '\u207f': 'n',
    }
    return ''.join(codepoints.get(char, char) for char in str(value))

def int_or_float(s):
    f = float(s)
    i = int(f)
    return i if i == f else f

class StripJunk(lark.Transformer):
    """
    Token stripper visitor class.

    This is done separately from the formula composer so that we can show the cleaned tree
    before debugging the conversion.
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
    def SUPERINT(self, token):
        """
        Return the integer value of a sequence of superscript digits.

        This is used in the charge rule as part of the valence specification for the atom.
        """
        return int(from_superscript(token.value))
    def DENSITYMODE(self, token):
        """
        Return the value of the DENSITYMODE token, either "n" or "i". If no mode is specified
        then a token value of None will be given to the density rule.
        """
        return token.value
    def CHARGE(self, token):
        """
        Return a sequence of plus and minus characters. By grammar rules they must all have
        the same sign.

        This is used in the charge rule as part of the valence specification for the atom.
        """
        return token.value
    def SUPERCHARGE(self, token):
        """
        Convert sequence of superscript plus and minus characters to ASCII plus and minus.

        This is used in the charge rule as part of the valence specification for the atom.
        """
        return from_superscript(token.value)
    def SYMBOL(self, token):
        """
        Look up the element in the periodic table and return it.

        Raise ValueError if the element doesn't exist.
        """
        try:
            return self._table.symbol(token.value)
        except Exception:
            raise ValueError(f"Element {token.value} doesn't exist")
    def FASTA(self, token):
        """
        Return the token value as the fasta sequence type: "dna", "rna" or "aa".
        """
        return token.value
    def SEQUENCE(self, token):
        """
        Return the token value as the fasta sequence string.
        """
        return token.value
    def fasta(self, tokens):
        """
        Return a fasta sequence and its type.

        Transform: [type, sequence] => ('fasta', type, sequence)
        """
        stype, sequence = tokens
        return 'fasta', stype, sequence
    def isotope(self, tokens):
        """
        Return the isotope number for the atom.

        Transform: [isotope] => isotope
        """
        return tokens[0]
    def charge(self, tokens):
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
        # Ca{++} => [None, '++'] = Ca.ion[2]
        # Ca{3--} => [3, '--'] = Ca.ion[-3]  # value has precedence over charge
        """
        print("in charge with", tokens)
        value, charge = tokens
        if value is None:
            value = len(charge)
        elif value and len(charge) > 1:
            self._raise_error(None, f"Using values of {value} for {value}{charge}")
        valence = value if charge[0] == '+' else -value
        return valence
    def atom(self, tokens):
        """
        Returns an atom from the periodic table.

        Usually this will use elements from the default table, but if an alternate table is
        provided to the ConvertTokens constructor then that will be used to retrieve the element
        from the symbol.

        Isotope and charge are optional. By using the rule "SYMBOL [isotope] [charge|supercharge]"
        with "[opt]" for the optional components rather "opt?", the missing components appear
        as None in the list of tokens. The "supercharge" option allows unicode superscripts to
        be used to specify charge rather than curly braces "{charge}".

        Raises an error if the symbol does not exist, does not have that isotope or doesn't
        allow that charge.

        Transform: ['symbol', isotope|None, charge|None] => atom

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

    def group(self, tokens):
        """
        Returns a sequence of (count, item) pairs, where item is an atom or a nested formula.
        Missing counts default to 1.

        Transform: [atom|formula, count|None, ...] => ((count, atom|formula), ...)
        """
        tokens = [1 if value is None else value for value in tokens]
        pairs = tuple((count, item) for item, count in zip(tokens[::2], tokens[1::2]))
        return pairs

    def composite(self, tokens):
        """
        Returns a sequence of (number, group) pairs. Each group is a sequence of (count, item)
        pairs, where item is an atom or a nested formula. Missing numbers default to 1.

        Transform: [number|None, group, ...] => ((number, group), ...) | ((count, atom), ...)

        Example CaCO3 6H2O: None, ((1, Ca), (1, C), (3, O)), 6, ((2, H), (1, O))]
        => ((1, ((1, Ca), (1, C), (3, O))), (6, ((2, H), (1, O))))

        Example CaCO3(H20)6: [[None, ((1, Ca), (1, C), (3, O), (6, formula('H2O')))]
        => ((1, Ca), (1, C), (3, O), (6, formula('H2O')))
        """
        # print("in composite", tokens)
        numbers = [1 if v is None else v for v in tokens[::2]]
        groups = tokens[1::2]
        pairs = tuple((number, group) for number, group in zip(numbers, groups))
        return pairs

    def fasta(self, tokens):
        """
        Returns the formula corresponding to the FASTA sequence, with the natural
        density set. Labile hydrogen use H[1] in the formula.

        The extra level of nesting in the return value is so that the fasta structure
        is like a composite with a single group containing a nested formula.

        Transform: [ /aa|dna|rna/, /[A-Z -*]+/ ] => (1, ((1, formula),))

        Example dna:CAGT: ['dna', 'CAGT'] => (1, ((1, C39H37H[1]10N15O25P4@1.69),))
        """
        # TODO: fasta is ignoring table when parsing
        # TODO: avoid circular imports
        # TODO: support other biochemicals (carbohydrate residues, lipids)
        from periodictable import fasta

        # print("in fasta", tokens)
        seq_type, seq = tokens
        if seq_type not in fasta.CODE_TABLES:
            raise ValueError(f"Invalid fasta sequence type '{seq_type}:'")
        seq = fasta.Sequence(name=None, sequence=seq, type=seq_type)
        group = ((1, seq.labile_formula),)
        composite = ((1, group),)
        return composite

    def density(self, tokens):
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

    def compound(self, tokens):
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

        Example NaCl@2.16i: [(1, ((1, Na), (1, Cl))), ('density', 2.16, 'i')] => NaCl@2.16i

        Example dna:CAGT: [((1, ((1, C39H37H[1]10N15O25P4@1.69n),)),), None] => C39H37H[1]10N15O25P4@1.69n

        Example CaCO3 6H2O: [((1, ((1, Ca), (1, C), (3, O))), (6, ((2, H), (1, O)))), None] => CaCO3(H2O)6

        Example CaCO3(H20)6: [((1, ((1, Ca), (1, C), (3, O), (6, H2O@None))),), None] => CaCO3(H2O)6
        """
        # print("in compound with", tokens)
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

        # print(f"compound = {formula} @ {formula.density}")
        return formula

    def weightpct(self, tokens):
        """
        Returns the percentage. The value has already be converted to a number.

        Used as the first percentage of a mix by weight mixture.

        Transform: [percent] => percent

        Example for "3 wt%": [3] => 3
        """
        return tokens[0]

    def volumepct(self, tokens):
        """
        Returns the percentage. The value has already be converted to a number.

        Used as the first percentage of a mix by volume mixture.

        Transform: [percent] => percent

        Example for "3 vol%": [3] => 3
        """
        return tokens[0]

    def percentage(self, tokens):
        """
        Returns the percentage. The value has already be converted to a number.

        Transform: [percent] => percent

        Example for " 3 % ": [3] => 3
        """
        return tokens[0]

    def byweight(self, tokens):
        """
        Returns mixture by wt% of the various components in the system.

        Raises ValueError if total exceeds 100%.

        Transform: [percent, formula, ..., percent, formula, formula] => formula

        Example: [76.95, D2O, H2O] => (D2O)3H2O
        """
        total = sum(tokens[:-1:2])
        if total > 100:
            raise ValueError(f"Total weight {total}% is more than 100%")
        pairs = [(compound, percent) for percent, compound in zip(tokens[:-1:2], tokens[1:-1:2])]
        pairs.append((tokens[-1], 100-total))
        # return 'byweight', [*pairs, last_pair]
        formula = _mix_by_weight_pairs(pairs)
        # print(f"byweight => {formula} @ {formula.density}")
        return formula

    def byvolume(self, tokens):
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
            raise ValueError(f"Total volume {total}% is more than 100%")
        pairs = [(compound, percent) for percent, compound in zip(tokens[:-1:2], tokens[1:-1:2])]
        pairs.append((tokens[-1], 100-total))
        # print("byvolume pairs", pairs)
        # print("byvolume density", [f.density for f, p in pairs])
        #return 'byvolume', pairs
        formula = _mix_by_volume_pairs(pairs)
        return formula

    def byamount(self, tokens):
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

    def layers(self, tokens):
        """
        Returns the mixture by volume of the various layers in the system.

        Raises ValueError if the density is missing from a component formula.

        Sets formula.thickness to the sum of the layer thicknesses.

        Transform: [quantity, formula, ...] => formula

        Example: [('length', 10.006, 'nm'), Ni, ('length', 3, 'mm'), Si] => NiSi164000
        """
        values = [value*LENGTH_UNITS[units] for dim, value, units in tokens[::2]]
        total = sum(values)
        percent = [(m/total)*100 for m in values]
        formula = _mix_by_volume_pairs(zip(tokens[1::2], percent))
        formula.thickness = total
        return formula

    def mixture(self, tokens):
        """
        Returns the formula representing the mixture, either byweight, byvolume, byamount or layers

        Transform: [formula] => formula
        """
        return tokens[0]

    def formula(self, tokens):
        """
        Return the formula representing the compound or mixture.

        Transform:  [formula] => formula
        """
        return tokens[0]

    def thickness(self, tokens):
        """
        Returns (dimension, value, unit) with dimension equal 'length'

        Transform: [value, ('length', unit)] => ('length', value, unit)

        Example: [5, ('length', 'nm')] => ('length', 5, 'nm')
        """
        value, (dim, units) = tokens
        return dim, value, units

    def quantity(self, tokens):
        """
        Returns (dimension, value, unit) with dimension equal 'mass' or 'volume'

        Transform: [value, (dimension, unit)] => (dimension, value, unit)

        Example: [5, ('mass', 'g')] => ('mass', 5, 'g')
        """
        value, (dim, units) = tokens
        return dim, value, units

    def start(self, tokens):
        """
        Return the final formula, with the original text attached.

        Sets formula.source to 'parse string' before returning.

        Transform: [formula] => formula
        """
        formula = tokens[0]
        # TODO: add the source string to the formula class attributes
        # Remember the string which was parsed
        formula.source = self._context
        return formula


def parse_formula(formula_str: str, table: PeriodicTable|None=None) -> Formula:
    """
    Parse a chemical formula, returning a structure with elements from the
    given periodic table.
    """
    cleanup = StripJunk()
    convert = ConvertTokens(formula_str, table=table)
    tree = parser(formula_str)
    tree = cleanup.transform(tree)
    tree = convert.transform(tree)
    return tree

examples = """
Co
dna:CAGT
(Co@5)
(((Co@5)@6))
CaCO3
CaCO₃
CaCO3+6H2O
CaCO3 6H2O
CaCO3(H2O)6
CaCO3 (H2O)6
(Ca(CO3)((H2O)6))
CaCO₃·6H₂O
DHO
!Ca{2++}  # could be interpreted as Ca{2+}
Ca⁺⁺  # also Ca{2+}
O²⁻
H[1]
H2O@1
D2O@1n
D2O @ 1.11
D2O@1.11i
HO{1-}
H[1]{1-}O
H2SO4
C3H4H[1]NO@1.29n
78.2H2O[16] + 21.8H2O[18] @1n
50 wt% Co // Ti
33 wt% Co // 33% Fe // Ti
! 93 wt% Co // 33% Fe // Ti  # More than 100%
20 vol% (10 wt% NaCl@2.16 // H2O@1) // D2O@1n
NaCl(H2O)29.1966(D2O)122.794@1.10i
5g NaCl // 50mL H2O@1
5g NaCl@2.16 // 50mL H2O@1
50 mL (45 mL H2O@1 // 5 g NaCl)@1.0707 // 20 mL D2O@1n
1 cm Si // 5 nm Cr // 10 nm Au
aa:RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV

# Error conditions. Mark with '!' so the exception is ignored
! Bl2Oh
! 5 Mg NaCl // 50mL H2O@1
! 4 nm NaCl@2.17// 50 g Si

"""

def check():
    cleanup = StripJunk()
    def filt(tree):
        #return tree
        tree = cleanup.transform(tree)
        # import pprint; pprint.pprint(tree)
        tree = convert.transform(tree)
        return tree

    for line in examples.split('\n'):
        formula = line.split('#')[0]
        bad = formula.startswith('!')
        if bad:
            formula = formula[1:]
        if formula:
            print(f"*** {line}")
            convert = ConvertTokens(text=formula)
            try:
                tree = filt(parser.parse(formula))
                #print(f" => {tree.pretty()}")
                density = getattr(tree, 'density', None)
                density_str = f" @ {density:.2f}" if density else ""
                print(f" => {tree}{density_str}")
                # TODO: structure not preserved in mixtures
                print(f"    {getattr(tree, 'structure', None)}")
            except Exception as exc:
                if bad:
                    print(f"!!! Error: {exc}")
                else:
                    raise

if __name__ == "__main__":
    check()