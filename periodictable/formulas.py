# This program is public domain
# Author: Paul Kienzle
"""
Chemical formula parser.
"""

from copy import copy
from math import pi, sqrt

# Requires that the pyparsing module is installed.

from pyparsing import (Literal, Optional, White, Regex,
                       ZeroOrMore, OneOrMore, Forward, StringEnd, Group)

from .core import default_table, isatom, isisotope, ision, change_table
from .constants import avogadro_number
from .util import cell_volume

PACKING_FACTORS = dict(cubic=pi/6, bcc=pi*sqrt(3)/8, hcp=pi/sqrt(18),
                       fcc=pi/sqrt(18), diamond=pi*sqrt(3)/16)


def mix_by_weight(*args, **kw):
    """
    Generate a mixture which apportions each formula by weight.

    :Parameters:
        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Returns:

        *formula* : Formula

    If density is not given, then it will be computed from the density
    of the components, assuming the components take up no more nor less
    space because they are in the mixture.  If component densities are
    not available, then the resulting density will not be computed.  The
    density calculation assumes the cell volume remains constant for the
    original materials, which is not in general the case.
    """
    table = default_table(kw.pop('table', None))
    density = kw.pop('density', None)
    natural_density = kw.pop('natural_density', None)
    name = kw.pop('name', None)
    if kw:
        raise TypeError("Unexpected arguments "+", ".join(kw.keys()))

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i], table=table), args[i+1])
             for i in range(0, len(args), 2)]
    result = _mix_by_weight_pairs(pairs)
    if natural_density:
        result.natural_density = natural_density
    if density:
        result.density = density
    if name:
        result.name = name
    return result

def _mix_by_weight_pairs(pairs):

    # Drop pairs with zero quantity
    # Note: must be first statement in order to accept iterators
    pairs = [(f, q) for f, q in pairs if q > 0]

    result = Formula()
    if pairs:
        # cell mass = mass
        # target mass = q
        # cell mass * n = target mass
        #   => n = target mass / cell mass
        #        = q / mass
        # scale this so that n = 1 for the smallest quantity
        scale = min(q/f.mass for f, q in pairs)
        for f, q in pairs:
            result += ((q/f.mass)/scale) * f
        if all(f.density for f, _ in pairs):
            volume = sum(q/f.density for f, q in pairs)/scale
            result.density = result.mass/volume
    return result

def mix_by_volume(*args, **kw):
    """
    Generate a mixture which apportions each formula by volume.

    :Parameters:
        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Returns:

        *formula* : Formula

    If density is not given, then it will be computed from the density
    of the components, assuming the components take up no more nor less
    space because they are in the mixture.  If component densities are
    not available, then a ValueError is raised. The  density calculation
    assumes the cell volume remains constant for the original materials,
    which is not in general the case.
    """
    table = default_table(kw.pop('table', None))
    density = kw.pop('density', None)
    natural_density = kw.pop('natural_density', None)
    name = kw.pop('name', None)
    if kw:
        raise TypeError("Unexpected arguments "+", ".join(kw.keys()))

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i], table=table), args[i+1])
             for i in range(0, len(args), 2)]
    result = _mix_by_volume_pairs(pairs)
    if natural_density:
        result.natural_density = natural_density
    if density:
        result.density = density
    if name:
        result.name = name
    return result

def _mix_by_volume_pairs(pairs):

    # Drop pairs with zero quantity
    # Note: must be first statement in order to accept iterators
    pairs = [(f, q) for f, q in pairs if q > 0]

    for f, _ in pairs:
        if f.density is None or f.density == 0.:
            raise ValueError("Need the mass density of "+str(f))

    result = Formula()
    if pairs:
        # cell volume = mass/density
        # target volume = q
        # cell volume * n = target volume
        #   => n = target volume / cell volume
        #        = q / (mass/density)
        #        = q * density / mass
        # scale this so that n = 1 for the smallest quantity
        scale = min(q*f.density/f.mass for f, q in pairs)
        for f, q in pairs:
            result += ((q*f.density/f.mass)/scale) * f

        volume = sum(q for _, q in pairs)/scale
        result.density = result.mass/volume

    return result

def formula(compound=None, density=None, natural_density=None,
            name=None, table=None):
    r"""
    Construct a chemical formula representation from a string, a
    dictionary of atoms or another formula.

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.

        *density* : float | |g/cm^3|
            Material density.  Not needed for single element formulas.

        *natural_density* : float | |g/cm^3|
            Material density assuming naturally occurring isotopes and no
            change in cell volume.

        *name* : string
            Common name for the molecule.

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Exceptions:
        *ValueError* : invalid formula initializer

    Example compounds:

       string:
          m = formula( "CaCO3+6H2O" )
       sequence of fragments:
          m = formula( [(1, Ca), (2, C), (3, O), (6, [(2, H), (1, O)]] )
       molecular math:
          m = formula( "CaCO3" ) + 6*formula( "H2O" )
       another formula (makes a copy):
          m = formula( formula("CaCO3+6H2O") )
       an atom:
          m = formula( Ca )
       nothing:
          m = formula()

    Operations:
       m.atoms returns {isotope: count, ...} for each atom in the compound.

    Formula strings consist of counts and atoms such as "CaCO3+6H2O".
    Groups can be separated by '+' or space, so "CaCO3 6H2O" works as well.
    Groups and be defined using parentheses, such as "CaCO3(H2O)6".
    Parentheses can nest: "(CaCO3(H2O)6)1"
    Isotopes are represented by index, e.g., "CaCO[18]3+6H2O".
    Counts can be integer or decimal, e.g. "CaCO3+(3HO0.5)2".
    Density can be specified in the formula using, e.g., "H2O@1". Isotopic
    formulas can use natural density, e.g., "D2O@1n", or the expected density
    with that isotope, e.g., "D2O@1.11".

    For full details see help(periodictable.formulas.formula_grammar)

    The chemical formula is designed for simple calculations such
    as molar mass, not for representing bonds or atom positions.
    However, we preserve the structure of the formula so that it can
    be used as a basis for a rich text representation such as
    matplotlib TeX markup.

    After creating a formula, a rough estimate of the density can be
    computed using::

        formula.density = formula.molecular_mass/formula.volume(packing_factor=...)

    The volume() calculation uses the covalent radii of the components and
    the known packing factor or crystal structure name.  If the lattice
    constants for the crystal are known, then they can be used instead::

        formula.density = formula.molecular_mass/formula.volume(a, b, c, alpha, beta, gamma)

    Formulas are designed for calculating quantities such as molar mass and
    scattering length density, not for representing bonds or atom positions.
    The representations are simple, but preserve some of the structure for
    display purposes.
    """
    if compound is None or compound == '':
        structure = tuple()
    elif isinstance(compound, Formula):
        structure = compound.structure
        if density is None and natural_density is None:
            density = compound.density
        if not name:
            name = compound.name
    elif isatom(compound):
        structure = ((1, compound), )
    elif isinstance(compound, dict):
        structure = _convert_to_hill_notation(compound)
    elif _is_string_like(compound):
        try:
            chem = parse_formula(compound, table=table)
            if name:
                chem.name = name
            if density is not None:
                chem.density = density
            elif natural_density is not None:
                chem.natural_density = natural_density
            return chem
        except ValueError as exception:
            raise ValueError(str(exception))
            #print "parsed", compound, "as", self
    else:
        try:
            structure = _immutable(compound)
        except:
            raise ValueError("not a valid chemical formula: "+str(compound))
    return Formula(structure=structure, name=name, density=density,
                   natural_density=natural_density)

class Formula:
    """
    Simple chemical formula representation.

    """
    def __init__(self, structure=tuple(), density=None, natural_density=None,
                 name=None):
        self.structure = structure
        self.name = name

        # If natural_density or density are specified, use them.
        # If only one element in the formula, use its density.
        # Otherwise, leave density unspecified, and let the user
        # assign it separately if they need it.
        if natural_density is not None:
            self.natural_density = natural_density
        elif density is not None:
            self.density = density
        elif len(self.atoms) == 1:
            # Note: density for isotopes already corrected for natural density
            atom = list(self.atoms.keys())[0]
            self.density = atom.density
        else:
            self.density = None

    @property
    def atoms(self):
        """
        { *atom*: *count*, ... }

        Composition of the molecule.  Referencing this attribute computes
        the *count* as the total number of each element or isotope in the
        chemical formula, summed across all subgroups.
        """
        return _count_atoms(self.structure)

    @property
    def hill(self):
        """
        Formula

        Convert the formula to a formula in Hill notation.  Carbon appears
        first followed by hydrogen then the remaining elements in alphabetical
        order.
        """
        return formula(self.atoms)

    def natural_mass_ratio(self):
        """
        Natural mass to isotope mass ratio.

        :Returns:
            *ratio* : float

        The ratio is computed from the sum of the masses of the individual
        elements using natural abundance divided by the sum of the masses
        of the isotopes used in the formula.  If the cell volume is
        preserved with isotope substitution, then the ratio of the masses
        will be the ratio of the densities.
        """
        total_natural_mass = total_isotope_mass = 0
        for el, count in self.atoms.items():
            try:
                natural_mass = el.element.mass
            except AttributeError:
                natural_mass = el.mass
            total_natural_mass += count * natural_mass
            total_isotope_mass += count * el.mass
        return total_natural_mass/total_isotope_mass

    @property
    def natural_density(self):
        """
        |g/cm^3|

        Density of the formula with specific isotopes of each element
        replaced by the naturally occurring abundance of the element
        without changing the cell volume.
        """
        return self.density*self.natural_mass_ratio()

    @natural_density.setter
    def natural_density(self, natural_density):
        self.density = natural_density / self.natural_mass_ratio()

    @property
    def mass(self):
        """
        atomic mass units u (C[12] = 12 u)

        Molar mass of the molecule.  Use molecular_mass to get the mass in
        grams.
        """
        mass = 0
        for el, count in self.atoms.items():
            mass += el.mass*count
        return mass

    @property
    def molecular_mass(self):
        """
        g

        Mass of the molecule in grams.
        """
        return self.mass/avogadro_number

    @property
    def charge(self):
        """
        Net charge of the molecule.
        """
        return sum([m*a.charge for a, m in self.atoms.items()])

    @property
    def mass_fraction(self):
        """
        Fractional mass representation of each element/isotope/ion.
        """
        total_mass = self.mass
        return dict((a, m*a.mass/total_mass) for a, m in self.atoms.items())

    def _pf(self):
        """
        packing factor  | unitless

        packing factor estimated from density.
        """
        return self.density

    def volume(self, *args, **kw):
        r"""
        Estimate unit cell volume.

        The crystal volume can be estimated from the element covalent radius
        and the atomic packing factor using::

            packing_factor = N_atoms V_atom / V_crystal

        Packing factors for a number of crystal lattice structures are defined.

        .. table:: Crystal lattice names and packing factors

         ======== ======================= ====================== ==============
         Code     Description             Formula                Packing factor
         ======== ======================= ====================== ==============
         cubic    simple cubic            $\pi/6$                0.52360
         bcc      body-centered cubic     $\pi\sqrt{3/8}$        0.68017
         hcp      hexagonal close-packed  $\pi/\sqrt{18}$        0.74048
         fcc      face-centered cubic     $\pi/\sqrt{18}$        0.74048
         diamond  diamond cubic           $\pi\sqrt{3/16}$       0.34009
         ======== ======================= ====================== ==============

        :Parameters:
            *packing_factor*  = 'hcp' : float or string
                Atomic packing factor.  If *packing_factor* is the name of
                a crystal lattice, use the *lattice* packing factor.
            *a*, *b*, *c* : float | |Ang|
                Lattice spacings. *b* and *c* default to *a*.
            *alpha*, *beta*, *gamma* : float | |deg|
                Lattice angles.  These default to 90\ |deg|

        :Returns:

            *volume* : float | |cm^3|
                Molecular volume.

        :Raises:

            *KeyError* : unknown lattice type

            *TypeError* : missing or bad lattice parameters

        Using the cell volume, mass density can be set with::

            formula.density = n*formula.molecular_mass/formula.volume()

        where n is the number of molecules per unit cell.

        Note: a single non-keyword argument is interpreted as a packing factor
        rather than a lattice spacing of 'a'.
        """
        # TODO: density estimated from H.covalent_radius is much too high
        #H_radius = kw.pop('H_radius', None)
        # Get packing factor
        if len(args) == 1 and not kw:
            packing_factor = args[0]
            args = []
        else:
            packing_factor = kw.pop('packing_factor', 'hcp')

        # Let cell_volume sort out its own parameters.
        if args or kw:
            return cell_volume(*args, **kw)*1e-24

        # Compute atomic volume
        V = 0
        for el, count in self.atoms.items():
            radius = el.covalent_radius
            #if el.number == 1 and H_radius is not None:
            #    radius = H_radius
            V += radius**3*count
        V *= 4.*pi/3

        # Translate packing factor from string
        try:
            _ = packing_factor + ""
        except Exception:
            pass
        else:
            packing_factor = PACKING_FACTORS[packing_factor.lower()]
        return V/packing_factor*1e-24

    def neutron_sld(self, *, wavelength=None, energy=None):
        """
        Neutron scattering information for the molecule.

        :Parameters:
            *wavelength* : float | |Ang|
                Wavelength of the neutron beam.

        :Returns:

            *sld* : (float, float, float) | |1e-6/Ang^2|
                Neutron scattering length density is returned as the tuple
                (*real*, *imaginary*, *incoherent*), or as (None, None, None)
                if the mass density is not known.

        .. deprecated:: 0.95
            Use periodictable.neutron_sld(formula) instead.
        """
        from .nsf import neutron_sld
        if self.density is None:
            return None, None, None
        return neutron_sld(self.atoms, density=self.density,
                           wavelength=wavelength, energy=energy)

    def xray_sld(self, *, energy=None, wavelength=None):
        """
        X-ray scattering length density for the molecule.

        :Parameters:
            *energy* : float | keV
                Energy of atom.

            *wavelength* : float | |Ang|
                Wavelength of atom.

            .. Note: One of *wavelength* or *energy* is required.

        :Returns:

            *sld* : (float, float) | |1e-6/Ang^2|
                X-ray scattering length density is returned as the tuple
                    (*real*, *imaginary*), or as (None, None) if the mass
                    density is not known.

        .. deprecated:: 0.95
            Use periodictable.xray_sld(formula) instead.
        """
        from .xsf import xray_sld
        if self.density is None:
            return None, None
        return xray_sld(self.atoms, density=self.density,
                        wavelength=wavelength, energy=energy)

    def change_table(self, table):
        """
        Replace the table used for the components of the formula.
        """
        self.structure = _change_table(self.structure, table)
        return self

    def replace(self, source, target, portion=1):
        """
        Create a new formula with one atom/isotope substituted for another.

        *formula* is the formula being updated.

        *source* is the isotope/element to be substituted.

        *target* is the replacement isotope/element.

        *portion* is the proportion of source which is substituted for target.
        """
        return _isotope_substitution(self, source, target, portion=portion)

    def __eq__(self, other):
        """
        Return True if two formulas represent the same structure. Note
        that they may still have different names and densities.
        Note: use hill representation for an order independent comparison.
        """
        if not isinstance(other, Formula):
            return False
        return self.structure == other.structure

    def __add__(self, other):
        """
        Join two formulas.
        """
        #print "adding", self, other
        if not isinstance(other, Formula):
            raise TypeError("expected formula+formula")
        ret = Formula()
        ret.structure = tuple(list(self.structure) + list(other.structure))
        return ret

    def __iadd__(self, other):
        """
        Extend a formula with another.
        """
        self.structure = tuple(list(self.structure) + list(other.structure))
        return self

    def __rmul__(self, other):
        """
        Provide a multiplier for formula.
        """
        #print "multiplying", self, other
        try:
            other += 0
        except TypeError:
            raise TypeError("n*formula expects numeric n")
        ret = copy(self)
        if other != 1 and self.structure:
            if len(self.structure) == 1:
                q, f = self.structure[0]
                ret.structure = ((other*q, f), )
            else:
                ret.structure = ((other, ret.structure), )
        return ret

    def __str__(self):
        return self.name if self.name else _str_atoms(self.structure)

    def __repr__(self):
        return "formula('%s')"%(str(self))


def _isotope_substitution(compound, source, target, portion=1):
    """
    Substitute one atom/isotope in a formula with another in some proportion.

    *compound* is the formula being updated.

    *source* is the isotope/element to be substituted.

    *target* is the replacement isotope/element.

    *portion* is the proportion of source which is substituted for target.
    """
    atoms = compound.atoms
    if source in atoms:
        mass = compound.mass
        mass_reduction = atoms[source]*portion*(source.mass - target.mass)
        density = compound.density * (mass - mass_reduction)/mass
        atoms[target] = atoms.get(target, 0) + atoms[source]*portion
        if portion == 1:
            del atoms[source]
        else:
            atoms[source] *= 1-portion
    else:
        density = compound.density
    return formula(atoms, density=density)


# TODO: Grammar should be independent of table
# TODO: Parser can't handle meters as 'm' because it conflicts with the milli prefix
LENGTH_UNITS = {'nm': 1e-9, 'um': 1e-6, 'mm': 1e-3, 'cm': 1e-2}
MASS_UNITS = {'ng': 1e-9, 'ug': 1e-6, 'mg': 1e-3, 'g': 1e+0, 'kg': 1e+3}
VOLUME_UNITS = {'nL': 1e-9, 'uL': 1e-6, 'mL': 1e-3, 'L': 1e+0}
LENGTH_RE = '('+'|'.join(LENGTH_UNITS.keys())+')'
MASS_VOLUME_RE = '('+'|'.join(list(MASS_UNITS.keys())+list(VOLUME_UNITS.keys()))+')'
def formula_grammar(table):
    """
    Construct a parser for molecular formulas.

    :Parameters:

        *table* = None : PeriodicTable
             If table is specified, then elements and their associated fields
             will be chosen from that periodic table rather than the default.

    :Returns:
        *parser* : pyparsing.ParserElement.
            The ``parser.parseString()`` method returns a list of
            pairs (*count, fragment*), where fragment is an *isotope*,
            an *element* or a list of pairs (*count, fragment*).

    """
    # TODO: fix circular imports
    # This ickiness is because the formula class returned from the circular
    # import of fasta does not match the local formula class.
    from .formulas import Formula

    # Recursive
    composite = Forward()
    mixture = Forward()

    # whitespace and separators
    space = Optional(White().suppress())
    separator = space+Literal('+').suppress()+space

    # Lookup the element in the element table
    symbol = Regex("[A-Z][a-z]?")
    symbol = symbol.setParseAction(lambda s, l, t: table.symbol(t[0]))

    # Translate isotope
    openiso = Literal('[').suppress()
    closeiso = Literal(']').suppress()
    isotope = Optional(~White()+openiso+Regex("[1-9][0-9]*")+closeiso,
                       default='0')
    isotope = isotope.setParseAction(lambda s, l, t: int(t[0]) if t[0] else 0)

    # Translate ion
    openion = Literal('{').suppress()
    closeion = Literal('}').suppress()
    ion = Optional(~White() +openion +Regex("([1-9][0-9]*)?[+-]") +closeion,
                   default='0+')
    ion = ion.setParseAction(lambda s, l, t: int(t[0][-1]+(t[0][:-1] if len(t[0]) > 1 else '1')))

    # Translate counts
    # TODO: regex should reject a bare '.' if we want to allow dots between formula parts
    fract = Regex("(0|[1-9][0-9]*|)([.][0-9]*)")
    fract = fract.setParseAction(lambda s, l, t: float(t[0]) if t[0] else 1)
    whole = Regex("(0|[1-9][0-9]*)")
    whole = whole.setParseAction(lambda s, l, t: int(t[0]) if t[0] else 1)
    number = Optional(~White()+(fract|whole), default=1)
    # TODO use unicode ₀₁₉ in the code below?
    sub_fract = Regex("(\u2080|[\u2081-\u2089][\u2080-\u2089]*|)([.][\u2080-\u2089]*)")
    sub_fract = sub_fract.setParseAction(lambda s, l, t: float(from_subscript(t[0])) if t[0] else 1)
    sub_whole = Regex("(\u2080|[\u2081-\u2089][\u2080-\u2089]*)")
    sub_whole = sub_whole.setParseAction(lambda s, l, t: int(from_subscript(t[0])) if t[0] else 1)
    sub_count = Optional(~White()+(fract|whole|sub_fract|sub_whole), default=1)

    # Fasta code
    fasta = Regex("aa|rna|dna") + Literal(":").suppress() + Regex("[A-Z *-]+")
    def convert_fasta(string, location, tokens):
        #print("fasta", string, location, tokens)
        # TODO: fasta is ignoring table when parsing
        # TODO: avoid circular imports
        # TODO: support other biochemicals (carbohydrate residues, lipids)
        from . import fasta
        seq_type, seq = tokens
        if seq_type not in fasta.CODE_TABLES:
            raise ValueError(f"Invalid fasta sequence type '{seq_type}:'")
        seq = fasta.Sequence(name=None, sequence=seq, type=seq_type)
        return seq.labile_formula
    fasta.setParseAction(convert_fasta)

    # Convert symbol, isotope, ion, count to (count, isotope)
    element = symbol+isotope+ion+sub_count
    def convert_element(string, location, tokens):
        """interpret string as element"""
        #print "convert_element received", tokens
        symbol, isotope, ion, count = tokens[0:4]
        if isotope != 0:
            symbol = symbol[isotope]
        if ion != 0:
            symbol = symbol.ion[ion]
        return (count, symbol)
    element = element.setParseAction(convert_element)

    # Convert "count elements" to a pair
    implicit_group = number+OneOrMore(element)
    def convert_implicit(string, location, tokens):
        """convert count followed by fragment"""
        #print "implicit", tokens
        count = tokens[0]
        fragment = tokens[1:]
        return fragment if count == 1 else (count, fragment)
    implicit_group = implicit_group.setParseAction(convert_implicit)

    # Convert "(composite) count" to a pair
    opengrp = space + Literal('(').suppress() + space
    closegrp = space + Literal(')').suppress() + space
    explicit_group = opengrp + composite + closegrp + sub_count
    def convert_explicit(string, location, tokens):
        """convert (fragment)count"""
        #print "explicit", tokens
        count = tokens[-1]
        fragment = tokens[:-1]
        return fragment if count == 1 else (count, fragment)
    explicit_group = explicit_group.setParseAction(convert_explicit)

    # Build composite from a set of groups
    group = implicit_group | explicit_group
    implicit_separator = separator | space
    composite << group + ZeroOrMore(implicit_separator + group)

    density = Literal('@').suppress() + number + Optional(Regex("[ni]"), default='i')
    compound = (composite|fasta) + Optional(density, default=None)
    def convert_compound(string, location, tokens):
        """convert material @ density or fasta @ density"""
        # Messiness: both composite and density can be one or more tokens
        # If density is missing then it is None, otherwise it is count + [ni]
        # Compound can be a sequence of (count, fragment) pairs, or if it is
        # a fasta sequence it may already be a formula.
        material = tokens[:-1] if tokens[-1] is None else tokens[:-2]
        #print("compound", material, type(material[0]), len(material))
        if len(material) == 1 and isinstance(material[0], Formula):
            formula = material[0]
        else:
            #print("unbundling material", material)
            formula = Formula(structure=_immutable(material))
        density, form = (None, None) if tokens[-1] is None else tokens[-2:]
        #if density is None and formula.density is None:
        #    # Estimate density from covalent radii and a 0.54 packing factor
        #    mass = formula.molecular_mass
        #    volume = formula.volume(packing_factor=0.54, H_radius=1.15)
        #    density, form = mass/volume, 'n'
        #    print(f"estimating density as {mass/volume=:.3f}")
        if form == 'n':
            formula.natural_density = density
        elif form == 'i':
            formula.density = density
        #print("compound", formula, f"{formula.density=:.3f}")
        return formula
    compound = compound.setParseAction(convert_compound)

    partsep = space + Literal('//').suppress() + space
    percent = Literal('%').suppress()
    weight = Regex("(w((eigh)?t)?|m(ass)?)").suppress()
    volume = Regex("v(ol(ume)?)?").suppress()
    weight_percent = (percent + weight) | (weight + percent) + space
    volume_percent = (percent + volume) | (volume + percent) + space
    by_weight = (number + weight_percent + mixture
                 + ZeroOrMore(partsep+number+(weight_percent|percent)+mixture)
                 + Optional(partsep + mixture, default=None))
    def _parts_by_weight_vol(tokens):
        #print("by weight or volume", tokens)
        if tokens[-1] is None:
            piece = tokens[1:-1:2]
            fract = [float(v) for v in tokens[:-1:2]]
            if abs(sum(fract) - 100) > 1e-12:
                raise ValueError(f"Formula percentages must sum to 100%, not {sum(fract)}")
        else:
            piece = tokens[1:-1:2] + [tokens[-1]]
            fract = [float(v) for v in tokens[:-1:2]]
            fract.append(100-sum(fract))
            if fract[-1] < 0:
                raise ValueError("Formula percentages must sum to less than 100%")
        #print piece, fract
        if len(piece) != len(fract):
            raise ValueError("Missing base component of mixture")
        return piece, fract
    def convert_by_weight(string, location, tokens):
        """convert mixture by wt% or mass%"""
        piece, fract = _parts_by_weight_vol(tokens)
        return _mix_by_weight_pairs(zip(piece, fract))
    mixture_by_weight = by_weight.setParseAction(convert_by_weight)

    by_volume = (number + volume_percent + mixture
                 + ZeroOrMore(partsep+number+(volume_percent|percent)+mixture)
                 + Optional(partsep + mixture, default=None))
    def convert_by_volume(string, location, tokens):
        """convert mixture by vol%"""
        piece, fract = _parts_by_weight_vol(tokens)
        return _mix_by_volume_pairs(zip(piece, fract))
    mixture_by_volume = by_volume.setParseAction(convert_by_volume)

    mixture_by_layer = Forward()
    layer_thick = Group(number + Regex(LENGTH_RE) + space)
    layer_part = (layer_thick + mixture) | (opengrp + mixture_by_layer + closegrp + sub_count)
    mixture_by_layer << layer_part + ZeroOrMore(partsep + layer_part)
    def convert_by_layer(string, location, tokens):
        """convert layer thickness '# nm material'"""
        if len(tokens) < 2:
            return tokens
        piece = []
        fract = []
        for p1, p2 in zip(tokens[0::2], tokens[1::2]):
            if isinstance(p1, Formula):
                f = p1.absthick * float(p2)
                p = p1
            else:
                f = float(p1[0]) * LENGTH_UNITS[p1[1]]
                p = p2
            piece.append(p)
            fract.append(f)
        total = sum(fract)
        vfract = [(v/total)*100 for v in fract]
        result = _mix_by_volume_pairs(zip(piece, vfract))
        result.thickness = total
        return result
    mixture_by_layer = mixture_by_layer.setParseAction(convert_by_layer)

    mixture_by_absmass = Forward()
    absmass_mass = Group(number + Regex(MASS_VOLUME_RE) + space)
    absmass_part = (absmass_mass + mixture) | (opengrp + mixture_by_absmass + closegrp + sub_count)
    mixture_by_absmass << absmass_part + ZeroOrMore(partsep + absmass_part)
    def convert_by_absmass(string, location, tokens):
        """convert mass '# mg material'"""
        if len(tokens) < 2:
            return tokens
        piece = []
        fract = []
        for p1, p2 in zip(tokens[0::2], tokens[1::2]):
            if isinstance(p1, Formula):
                p = p1
                f = p1.total_mass * float(p2)
            else:
                p = p2
                value = float(p1[0])
                if p1[1] in VOLUME_UNITS:
                    # convert to volume in liters to mass in grams before mixing
                    if p.density is None:
                        raise ValueError("Need the mass density of "+str(p))
                    f = value * VOLUME_UNITS[p1[1]] * 1000.*p.density
                else:
                    f = value * MASS_UNITS[p1[1]]
            piece.append(p)
            fract.append(f)

        total = sum(fract)
        mfract = [(m/total)*100 for m in fract]
        result = _mix_by_weight_pairs(zip(piece, mfract))
        result.total_mass = total
        return result
    mixture_by_absmass = mixture_by_absmass.setParseAction(convert_by_absmass)

    ungrouped_mixture = (mixture_by_weight | mixture_by_volume
                         | mixture_by_layer | mixture_by_absmass)
    grouped_mixture = opengrp + ungrouped_mixture + closegrp + Optional(density, default=None)
    def convert_mixture(string, location, tokens):
        """convert (mixture) @ density"""
        formula = tokens[0]
        if tokens[-1] == 'n':
            formula.natural_density = tokens[-2]
        elif tokens[-1] == 'i':
            formula.density = tokens[-2]
        # elif tokens[-1] is None
        return formula
    grouped_mixture = grouped_mixture.setParseAction(convert_mixture)

    mixture << (compound | grouped_mixture)
    formula = (compound | ungrouped_mixture | grouped_mixture)
    grammar = Optional(formula, default=Formula()) + StringEnd()

    grammar.setName('Chemical Formula')
    return grammar

_PARSER_CACHE = {}
def parse_formula(formula_str, table=None):
    """
    Parse a chemical formula, returning a structure with elements from the
    given periodic table.
    """
    table = default_table(table)
    if table not in _PARSER_CACHE:
        _PARSER_CACHE[table] = formula_grammar(table)
    parser = _PARSER_CACHE[table]
    #print(parser)
    return parser.parseString(formula_str)[0]

def _count_atoms(seq):
    """
    Traverse formula structure, counting the total number of atoms.
    """
    total = {}
    for count, fragment in seq:
        if isinstance(fragment, (list, tuple)):
            partial = _count_atoms(fragment)
        else:
            partial = {fragment: 1}
        for atom, atom_count in partial.items():
            if atom not in total:
                total[atom] = 0
            total[atom] += atom_count*count
    return total

def count_elements(compound, by_isotope=False):
    """
    Element composition of the molecule.

    Returns {*element*: *count*, ...} where the *count* is the total number
    of each element in the chemical formula, summed across all isotopes and
    ionization levels.

    If *by_isotope* is True, then sum across ionization
    levels, keeping the individual isotopes separate.
    """
    total = {}
    # Note: could accumulate charge at the same time as counting elements.
    for part, count in formula(compound).atoms.items():
        # Resolve isotopes and ions to the underlying element. Four cases:
        #    isotope with charge needs fragment.element.element
        #    isotope without charge needs fragment.element
        #    element with charge needs fragment.element
        #    element without charge needs fragment
        if ision(part):
            part = part.element
        if not by_isotope:
            part = getattr(part, "element", part)
        total[part] = count + total.get(part, 0)
    return total

def _immutable(seq):
    """
    Traverse formula structure, checking that the counts are numeric and
    units are atoms.  Returns an immutable copy of the structure, with all
    lists replaced by tuples.
    """
    if isatom(seq):
        return seq
    return tuple((count+0, _immutable(fragment)) for count, fragment in seq)

def _change_table(seq, table):
    """Converts lists to tuples so that structure is immutable."""
    if isatom(seq):
        return change_table(seq, table)
    return tuple((count, _change_table(fragment, table))
                 for count, fragment in seq)

def _hill_compare(a, b):
    """
    Compare elements in standard order.
    """
    if a.symbol == b.symbol:
        a = a.isotope if isisotope(a) else 0
        b = b.isotope if isisotope(b) else 0
        return cmp(a, b)
    elif a.symbol in ("C", "H"):
        if b.symbol in ("C", "H"):
            return cmp(a.symbol, b.symbol)
        else:
            return -1
    else:
        if b.symbol in ("C", "H"):
            return 1
        else:
            return cmp(a.symbol, b.symbol)

def _hill_key(a):
    return "".join((("0" if a.symbol in ("C", "H") else "1"),
                    a.symbol,
                    "%4d"%(a.isotope if isisotope(a) else 0)))

def _convert_to_hill_notation(atoms):
    """
    Return elements listed in standard order.
    """
    #return [(atoms[el], el) for el in sorted(atoms.keys(), cmp=_hill_compare)]
    return [(atoms[el], el) for el in sorted(atoms.keys(), key=_hill_key)]

def _str_one_atom(fragment):
    # Normal isotope string form is #-Yy, but we want Yy[#]
    if isisotope(fragment) and 'symbol' not in fragment.__dict__:
        ret = "%s[%d]"%(fragment.symbol, fragment.isotope)
    else:
        ret = fragment.symbol
    if fragment.charge != 0:
        sign = '+' if fragment.charge > 0 else '-'
        value = str(abs(fragment.charge)) if abs(fragment.charge) > 1 else ''
        ret += '{'+value+sign+'}'
    return ret

def _str_atoms(seq):
    """
    Convert formula structure to string.
    """
    #print "str", seq
    ret = ""
    for count, fragment in seq:
        if isatom(fragment):
            ret += _str_one_atom(fragment)
            if count != 1:
                ret += "%g"%count
        else:
            if count == 1:
                piece = _str_atoms(fragment)
            else:
                piece = "(%s)%g"%(_str_atoms(fragment), count)
            #ret = ret+" "+piece if ret else piece
            ret += piece

    return ret

def _is_string_like(val):
    """Returns True if val acts like a string"""
    try:
        val+''
    except Exception:
        return False
    return True

def from_subscript(value):
    subscript_codepoints = {
        '\u2080': '0', '\u2081': '1', '\u2082': '2', '\u2083': '3',
        '\u2084': '4', '\u2085': '5', '\u2086': '6', '\u2087': '7',
        '\u2088': '8', '\u2089': '9', '\u208a': '+', '\u208b': '-',
        '\u208c': '=', '\u208d': '(', '\u208e': ')',

        '\u2090': 'a', '\u2091': 'e', '\u2092': 'o', '\u2093': 'x',
        '\u2095': 'h', '\u2096': 'k', '\u2097': 'l',
        '\u2098': 'm', '\u2099': 'n', '\u209a': 'p', '\u209b': 's',
        '\u209c': 't',
    }
    return ''.join(subscript_codepoints.get(char, char) for char in str(value))

def unicode_subscript(value):
    # Unicode subscript codepoints. Note that decimal point looks okay as subscript
    subscript_codepoints = {
        '0': '\u2080', '1': '\u2081', '2': '\u2082', '3': '\u2083',
        '4': '\u2084', '5': '\u2085', '6': '\u2086', '7': '\u2087',
        '8': '\u2088', '9': '\u2089', '+': '\u208a', '-': '\u208b',
        '=': '\u208c', '(': '\u208d', ')': '\u208e',

        'a': '\u2090', 'e': '\u2091', 'o': '\u2092', 'x': '\u2093',
        'h': '\u2095', 'k': '\u2096', 'l': '\u2097',
        'm': '\u2098', 'n': '\u2099', 'p': '\u209a', 's': '\u209b',
        't': '\u209c',

        '\u2013': '\u208b', # en-dash is same as dash
        '\u2014': '\u208b', # em-dash is same as dash
    }
    return ''.join(subscript_codepoints.get(char, char) for char in str(value))

def unicode_superscript(value):
    # Unicode subscript codepoints. Note that decimal point looks okay as subscript
    superscript_codepoints = {
        #'.': '\u00B0',  # degree symbol looks too much like zero
        #'.': ' \u02D9',  # dot above modifier looks okay in a floating string, but risky
        #'.': ' \u0307',  # space with dot above?
        #'.': '\u22C5', # math dot operator
        '.': '\u1427',  # Canadian aboriginal extended block dot (looks good on mac)
        '2': '\u00B2', '3': '\u00B3',
        '1': '\u00B9',
        '0': '\u2070', 'i': '\u2071',
        '4': '\u2074', '5': '\u2075', '6': '\u2076', '7': '\u2077',
        '9': '\u2078', '0': '\u2079', '+': '\u207a', '-': '\u207b',
        '=': '\u207c', '(': '\u207d', ')': '\u207e', 'n': '\u207f',

        '\u2013': '\u207b', # en-dash is same as dash
        '\u2014': '\u207b', # em-dash is same as dash
    }
    return ''.join(superscript_codepoints.get(char, char) for char in str(value))

SUBSCRIPT = {
    # The latex renderer should work for github style markdown
    'latex': lambda text: f'$_{{{text}}}$',
    'html': lambda text: f'<sub>{text}</sub>',
    'unicode': unicode_subscript,
    'plain': lambda text: text
}
def pretty(compound, mode='unicode'):
    """
    Convert the formula to a string. The *mode* can be 'unicode', 'html' or
    'latex' depending on how subscripts should be rendered. If *mode* is 'plain'
    then don't use subscripts for the element quantities.

    Use *pretty(compound.hill)* for a more compact representation.
    """
    return _pretty(compound.structure, SUBSCRIPT[mode])

def _pretty(structure, subscript):
    # TODO: if superscript is not None then render O[16] as {}^{16}O
    parts = []
    for count, part in structure:
        if isinstance(part, tuple):
            if count == 1:
                parts.append(_pretty(part, subscript))
            else:
                parts.append(f'({_pretty(part, subscript)}){subscript(count)}')
        elif count == 1:
            parts.append(f'{_str_one_atom(part)}')
        else:
            parts.append(f'{_str_one_atom(part)}{subscript(count)}')
    return ''.join(parts)


def demo():
    import sys
    compound = formula(sys.argv[1])
    if compound.density is None:
        print(f"Missing density for {pretty(compound.hill)}")
    else:
        print(f"{pretty(compound.hill)}@{compound.density:.3f} sld={compound.neutron_sld()}")
        #print(pretty(compound.hill, 'latex'))
        #print(pretty(compound.hill, 'html'))
        #print(pretty(compound.hill, 'unicode'))
        #print(pretty(compound.hill, 'plain'))

if __name__ == "__main__":
    demo()
