# This program is public domain
# Author: Paul Kienzle
"""
Core classes for the periodic table.

* :class:`PeriodicTable`
   The periodic table with attributes for each element.

   .. Note:: PeriodicTable is not a singleton class.  Use ``periodictable.element``
      to access the common table.

* :class:`Element`
   Element properties such as name, symbol, mass, density, etc.

* :class:`Isotope`
   Isotope properties such as mass, density and neutron scattering factors.

* :class:`Ion`
   Ion properties such as charge.

Elements are accessed from a periodic table using ``table[number]``,
``table.name`` or ``table.symbol`` where *symbol* is the two letter symbol.
Individual isotopes are accessed using ``element[isotope]``. Individual ions
are references using ``element.ion[charge]``.  Note that
``element[isotope].ion[charge].mass`` will depend on the particular charge
since we subtract the charge times the rest mass of the electron from the
overall mass.

Helper functions:

* :func:`delayed_load`
    Delay loading the element attributes until they are needed.

* :func:`get_data_path`
    Return the path to the periodic table data files.

* :func:`define_elements`
    Define external variables for each element in namespace.

* :func:`isatom`, :func:`iselement`, :func:`isisotope`, :func:`ision`
    Tests for different types of structure components.

* :func:`default_table`
    Returns the common periodic table.

* :func:`change_table`
    Return the same item from a different table.

.. seealso::

    :ref:`Adding properties <extending>` for details on extending the periodic
    table with your own attributes.

    :ref:`Custom tables <custom-table>` for details on managing your own
    periodic table with custom values for the attributes.

"""

__docformat__ = 'restructuredtext en'
__all__ = ['delayed_load', 'define_elements', 'get_data_path',
           'default_table', 'change_table',
           'Ion', 'Isotope', 'Element', 'PeriodicTable',
           'isatom', 'iselement', 'isisotope', 'ision']

from pathlib import Path
from typing import TYPE_CHECKING, cast, Any, Union, TypeVar, List, overload, Literal
from collections.abc import Sequence, Callable, Iterator

if TYPE_CHECKING:
    from .nsf import Neutron
    from .activation import ActivationResult
    from .xsf import Xray
    from .magnetic_ff import MagneticFormFactor

from . import constants

Atom = Union["Element", "Isotope", "Ion"]
AtomVar = TypeVar('AtomVar', bound=Atom)

PUBLIC_TABLE_NAME = "public"

def delayed_load(all_props: Sequence[str], loader: Callable[[], None], element=True, isotope=False, ion=False):
    """
    Delayed loading of an element property table.  When any of property
    is first accessed the loader will be called to load the associated
    data. The help string starts out as the help string for the loader
    function. The attribute may be associated with any of :class:`Isotope`,
    :class:`Ion`, or :class:`Element`. Some properties, such as
    :mod:`mass <periodictable.mass>`, have both an isotope property for the
    mass of specific isotopes, as well as an element property for the
    mass of the collection of isotopes at natural abundance.  Set the
    keyword flags *element*, *isotope* and/or *ion* to specify which
    of these classes will be assigned specific information on load.
    """
    def clearprops():
        """
        Remove the properties so that the attribute can be accessed
        directly.
        """
        if element:
            for p in all_props:
                delattr(Element, p)
        if isotope:
            for p in all_props:
                delattr(Isotope, p)
        if ion:
            for p in all_props:
                delattr(Ion, p)

    def getter(propname):
        """
        Property getter for attribute propname.

        The first time the prop is accessed, the prop itself will be
        deleted and the data loader for the property will be called
        to set the real values.  Subsequent references to the property
        will be to the actual data.
        """
        def getfn(el):
            #print "get", el, propname
            clearprops()
            loader()
            return getattr(el, propname)
        return getfn

    def setter(propname):
        """
        Property setter for attribute propname.

        This function is assumed to be called when the data loader for the
        attribute is called before the property is referenced (for example,
        if somebody imports periodictable.xsf before referencing Ni.xray).
        In this case, we simply need to clear the delayed load property and
        let the loader set the values as usual.

        If the user tries to override a value in the table before first
        referencing the table, then the above assumption is false. E.g.,
        "Ni.K_alpha=5" followed by "print Cu.K_alpha" will yield an
        undefined Cu.K_alpha. This will be difficult for future users
        to debug.
        """
        def setfn(el, value):
            #print "set", el, propname, value
            clearprops()
            # Since the property is now cleared the lazy loader will not be
            # triggered in the getter for a different element, so call it here.
            loader()
            setattr(el, propname, value)
        return setfn

    doc = loader.__doc__
    if element:
        for p in all_props:
            prop = property(getter(p), setter(p), doc=doc)
            setattr(Element, p, prop)

    if isotope:
        for p in all_props:
            prop = property(getter(p), setter(p), doc=doc)
            setattr(Isotope, p, prop)

    if ion:
        for p in all_props:
            prop = property(getter(p), setter(p), doc=doc)
            setattr(Ion, p, prop)

# Define the element names from the element table.
class PeriodicTable:
    """
    Defines the periodic table of the elements with isotopes.
    Individidual elements are accessed by name, symbol or atomic number.
    Individual isotopes are addressable by ``element[mass_number]`` or
    ``elements.isotope(element name)``, ``elements.isotope(element symbol)``.

    For example, the following all retrieve iron:

    .. doctest::

        >>> from periodictable import *
        >>> print(elements[26])
        Fe
        >>> print(elements.Fe)
        Fe
        >>> print(elements.symbol('Fe'))
        Fe
        >>> print(elements.name('iron'))
        Fe
        >>> print(elements.isotope('Fe'))
        Fe


    To get iron-56, use:

    .. doctest::

        >>> print(elements[26][56])
        56-Fe
        >>> print(elements.Fe[56])
        56-Fe
        >>> print(elements.isotope('56-Fe'))
        56-Fe


    Deuterium and tritium are defined as 'D' and 'T'.

    To show all the elements in the table, use the iterator:

    .. doctest::

        >>> from periodictable import *
        >>> for el in elements:  # lists the element symbols
        ...     print("%s %s"%(el.symbol, el.name))  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        H hydrogen
        He helium
        ...
        Og oganesson


    .. Note::
           Properties can be added to the elements as needed, including *mass*,
           *nuclear* and *X-ray* scattering cross sections.
           See section :ref:`Adding properties <extending>` for details.
    """
    properties: List[str]  # list method shadows builtin list, so using List instead
    """Properties loaded into the table"""

    # Tedious listing of available elements for typed table.El access
    n: "Element"
    H: "Element"
    D: "Isotope"
    T: "Isotope"
    He: "Element"
    Li: "Element"
    Be: "Element"
    B: "Element"
    C: "Element"
    N: "Element"
    O: "Element"
    F: "Element"
    Ne: "Element"
    Na: "Element"
    Mg: "Element"
    Al: "Element"
    Si: "Element"
    P: "Element"
    S: "Element"
    Cl: "Element"
    Ar: "Element"
    K: "Element"
    Ca: "Element"
    Sc: "Element"
    Ti: "Element"
    V: "Element"
    Cr: "Element"
    Mn: "Element"
    Fe: "Element"
    Co: "Element"
    Ni: "Element"
    Cu: "Element"
    Zn: "Element"
    Ga: "Element"
    Ge: "Element"
    As: "Element"
    Se: "Element"
    Br: "Element"
    Kr: "Element"
    Rb: "Element"
    Sr: "Element"
    Y: "Element"
    Zr: "Element"
    Nb: "Element"
    Mo: "Element"
    Tc: "Element"
    Ru: "Element"
    Rh: "Element"
    Pd: "Element"
    Ag: "Element"
    Cd: "Element"
    In: "Element"
    Sn: "Element"
    Sb: "Element"
    Te: "Element"
    I: "Element"
    Xe: "Element"
    Cs: "Element"
    Ba: "Element"
    La: "Element"
    Ce: "Element"
    Pr: "Element"
    Nd: "Element"
    Pm: "Element"
    Sm: "Element"
    Eu: "Element"
    Gd: "Element"
    Tb: "Element"
    Dy: "Element"
    Ho: "Element"
    Er: "Element"
    Tm: "Element"
    Yb: "Element"
    Lu: "Element"
    Hf: "Element"
    Ta: "Element"
    W: "Element"
    Re: "Element"
    Os: "Element"
    Ir: "Element"
    Pt: "Element"
    Au: "Element"
    Hg: "Element"
    Tl: "Element"
    Pb: "Element"
    Bi: "Element"
    Po: "Element"
    At: "Element"
    Rn: "Element"
    Fr: "Element"
    Ra: "Element"
    Ac: "Element"
    Th: "Element"
    Pa: "Element"
    U: "Element"
    Np: "Element"
    Pu: "Element"
    Am: "Element"
    Cm: "Element"
    Bk: "Element"
    Cf: "Element"
    Es: "Element"
    Fm: "Element"
    Md: "Element"
    No: "Element"
    Lr: "Element"
    Rf: "Element"
    Db: "Element"
    Sg: "Element"
    Bh: "Element"
    Hs: "Element"
    Mt: "Element"
    Ds: "Element"
    Rg: "Element"
    Cn: "Element"
    Nh: "Element"
    Fl: "Element"
    Mc: "Element"
    Lv: "Element"
    Ts: "Element"
    Og: "Element"

    def __init__(self, table: str) -> None:
        if table in PRIVATE_TABLES:
            raise ValueError("Periodic table '%s' is already defined"%table)
        PRIVATE_TABLES[table] = self
        self.properties = []
        self._element = {}
        for Z, (name, symbol, ions, uncommon_ions) in element_base.items():
            element = Element(name=name.lower(), symbol=symbol, Z=Z,
                              ions=tuple(sorted(ions+uncommon_ions)), table=table)
            self._element[element.number] = element
            setattr(self, symbol, element)
            PeriodicTable.__annotations__[symbol] = Element

        # There are two specially named isotopes D and T
        self.D = self.H.add_isotope(2)
        self.D.name = 'deuterium'
        self.D.symbol = 'D'
        self.T = self.H.add_isotope(3)
        self.T.name = 'tritium'
        self.T.symbol = 'T'

    def __getitem__(self, Z: int) -> "Element":
        """
        Retrieve element Z.
        """
        return self._element[Z]

    def __iter__(self) -> Iterator["Element"]:
        """
        Process the elements in Z order
        """
        # CRUFT: Since 3.7 dictionaries use insertion order, so no need to sort
        elements = sorted(self._element.items())
        # Skipping the first entry (neutron) in the iterator
        for _, el in elements[1:]:
            yield el

    @overload
    def symbol(self, input: Literal["D"]) -> "Isotope": ... # type: ignore[overload-overlap]

    @overload
    def symbol(self, input: Literal["T"]) -> "Isotope": ... # type: ignore[overload-overlap]

    @overload
    def symbol(self, input: str) -> "Element": ...

    def symbol(self, input: str) -> Union["Element", "Isotope"]:
        """
        Lookup the an element in the periodic table using its symbol.  Symbols
        are included for 'D' and 'T', deuterium and tritium.

        :Parameters:
            *input* : string
                Element symbol to be looked up in periodictable.

        :Returns: Element

        :Raises:
            ValueError if the element symbol is not defined.

        For example, print the element corresponding to 'Fe':

        .. doctest::

            >>> import periodictable
            >>> print(periodictable.elements.symbol('Fe'))
            Fe
        """
        if hasattr(self, input):
            value = getattr(self, input)
            if isinstance(value, (Element, Isotope)):
                return value
        raise ValueError("unknown element "+input)

    def name(self, input: str) -> Union["Element", "Isotope"]:
        """
        Lookup an element given its name.

        :Parameters:
            *input* : string
                Element name to be looked up in periodictable.

        :Returns: Element (or Isotope for "D" and "T")

        :Raises:
            *ValueError* if element does not exist.

        For example, print the element corresponding to 'iron':

        .. doctest::

            >>> import periodictable
            >>> print(periodictable.elements.name('iron'))
            Fe
        """
        for el in self:
            if input == el.name:
                return el
        if input == self.D.name:
            return self.D
        if input == self.T.name:
            return self.T
        raise ValueError("unknown element "+input)

    def isotope(self, input: str) -> Union["Element", "Isotope"]:
        """
        Lookup the element or isotope in the periodic table. Elements
        are assumed to be given by the standard element symbols. Isotopes
        are given by number-symbol, or 'D' and 'T' for 2-H and 3-H.

        :Parameters:
            *input* : string
                Element name or isotope to be looked up in periodictable.

        :Returns: Element

        :Raises:
            *ValueError* if element or isotope is not defined.

        For example, print the element corresponding to '58-Ni'.

        .. doctest::

            >>> import periodictable
            >>> print(periodictable.elements.isotope('58-Ni'))
            58-Ni
        """
        # Parse #-Sym or Sym
        # If # is not an integer, set isotope to -1 so that the isotope
        # lookup will fail later.
        parts = input.split('-')
        if len(parts) == 1:
            isotope = 0
            symbol = parts[0]
        elif len(parts) == 2:
            try:
                isotope = int(parts[0])
            except Exception:
                isotope = -1
            symbol = parts[1]
        else:
            symbol = ''
            isotope = -1

        # All elements are attributes of the table
        # Check that the attribute is an Element or an Isotope (for D or T)
        # If it is an element, check that the isotope exists
        if hasattr(self, symbol):
            attr = getattr(self, symbol)
            if isinstance(attr, Element):
                # If no isotope, return the element
                if isotope == 0:
                    return attr
                # If isotope, check that it is valid
                if isotope in attr.isotopes:
                    return attr[isotope]
            elif isinstance(attr, Isotope):
                # D, T must not have an associated isotope; 4-D is meaningless.
                if isotope == 0:
                    return attr

        # If we can't parse the string as an element or isotope, raise an error
        raise ValueError("unknown element "+input)

    # TODO: list method shadows builtin list in "properties: list[str]" above
    def list(self, *props, **kw) -> None:
        """
        Print a list of elements with the given set of properties.

        :Parameters:
            *prop1*, *prop2*, ... : string
                Name of the properties to print
            *format*: string
                Template for displaying the element properties, with one
                % for each property.

        :Returns: None

        For example, print a table of mass and density.

        .. doctest::

            >>> from periodictable import elements
            >>> elements.list('symbol', 'mass', 'density',
            ...     format="%-2s: %6.2f u %6.2f g/cm^3") # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
            H :   1.01 u   0.07 g/cm^3
            He:   4.00 u   0.12 g/cm^3
            Li:   6.94 u   0.53 g/cm^3
            ...
            Bk: 247.00 u  14.00 g/cm^3
        """
        #TODO: accept template strings
        #TODO: override signature in sphinx with
        #    .. method:: list(prop1, prop2, ..., format='')
        format = kw.pop('format', None)
        assert not kw  # make sure *format* is the only keyword argument
        for el in self:
            try:
                L = tuple(getattr(el, p) for p in props)
            except AttributeError:
                # Skip elements which don't define all the attributes
                continue
            # Skip elements with a value of None
            if any(v is None for v in L):
                continue

            if format is None:
                print(" ".join(str(p) for p in L))
            else:
                #try:
                print(format%L)
                #except:
                #    print "format", format, "args", L
                #    raise

# TODO: types for properties are not being handled correctly
# I've left them as simple "name: float" for now so that the editor will see them.
# Because I'm using delegation via __getattr__ the properties need to be listed in
# _AtomBase, but I can't define them as properties there because it confuses the
# delegation mechanism. Possible work-arounds are to build the delegation into the
# property, or figure out how to define a property return type that is seen by the
# static type analyzers without defining the property itself.
class _AtomBase:
    """
    Attributes common to element, isotope and ion.

    This class is defined only for type hinting. Some of the attributes are accessible
    as properties. Those not defined in isotope or ion are delegated to the base element
    through attribute access magic.
    """
    # attributes delegated to Element class
    table: str
    name: str
    symbol: str
    number: int
    ions: tuple[int, ...]
    ion: "IonSet" # TODO: could be IonSet["Element"] or IonSet["Isotope"]
    charge: int # element (=0), isotope (delegate to element), ion (!= 0)

    #element: Union["Element", "Isotope"] # ion or isotope

    # mass.py
    mass: float # property
    _mass: float # internal
    _mass_unc: float # Not yet official, but the data is loaded for some tables
    mass_units: str # element, isotope, ion
    #abundance: float  # isotope only
    #abundance_units: str # isotope only

    # covalent_radius.py
    covalent_radius: float|None # element
    covalent_radius_uncertainty: float|None # element
    covalent_radius_units: str # element

    # crystal_structure.py
    crystal_structure: dict[str, Any]|None # element

    # density.py
    density: float # property
    _density: float
    density_units: str
    interatomic_distance: float # property
    interatomic_distance_units: str
    number_density: float # property
    number_density_units: str
    density_caveat: str|None

    # nsf.py
    neutron: "Neutron" # element and isotope
    #nuclear_spin: str # isotope only

    # activation.py
    #neutron_activation: tuple["ActivationResult"]|None  # isotope only

    # xsf.py
    K_alpha: float|None # element
    K_beta1: float|None # element
    K_alpha_units: str # element
    K_beta1_units: str # element
    xray: "Xray" # element

    # magnetic_ff.py
    magnetic_ff: dict[int, "MagneticFormFactor"] # element

class Ion(_AtomBase):
    """
    Periodic table entry for an individual ion.

    An ion is associated with an element. In addition to the element
    properties (*symbol*, *name*, *atomic number*), it has specific ion
    properties (*charge*). Properties not specific to the ion (i.e., *charge*)
    are retrieved from the associated element.
    """
    element: Union["Element", "Isotope"]
    # charge: int # inherited from _AtomBase

    # TODO: abundance and activation need to be defined for charged isotopes.

    def __init__(self, element: Union["Element", "Isotope"], charge: int):
        self.element = element
        self.charge = charge
    def __getattr__(self, attr: str) -> Any:
        return getattr(self.element, attr)
    @property
    def mass(self) -> float: # type: ignore[override]
        return getattr(self.element, 'mass') - constants.electron_mass*self.charge
    def __str__(self) -> str:
        sign = '+' if self.charge > 0 else '-'
        value = '%d'%abs(self.charge) if abs(self.charge) > 1 else ''
        charge_str = '{'+value+sign+'}' if self.charge != 0 else ''
        return str(self.element)+charge_str
    def __repr__(self) -> str:
        return repr(self.element)+'.ion[%d]'%self.charge
    def __reduce__(self):
        if isinstance(self.element, Isotope):
            return _make_isotope_ion, (self.element.table,
                                       self.element.number,
                                       self.element.isotope,
                                       self.charge)
        else:
            return _make_ion, (self.element.table,
                               self.element.number,
                               self.charge)

class IonSet:
    element_or_isotope: Union["Element", "Isotope"]
    ionset: dict[int, "Ion"]

    def __init__(self, element_or_isotope: Union["Element", "Isotope"]):
        self.element_or_isotope = element_or_isotope
        self.ionset = {}

    def __getitem__(self, charge: int) -> Ion:
        if charge not in self.ionset:
            if charge not in self.element_or_isotope.ions:
                raise ValueError("%(charge)d is not a valid charge for %(symbol)s"
                                 % dict(charge=charge,
                                        symbol=self.element_or_isotope.symbol))
            self.ionset[charge] = Ion(self.element_or_isotope, charge)
        return self.ionset[charge]

class Isotope(_AtomBase):
    """
    Periodic table entry for an individual isotope.

    An isotope is associated with an element.  In addition to the element
    properties (*symbol*, *name*, *atomic number*), it has specific isotope
    properties (*isotope number*, *nuclear spin*, *relative abundance*).
    Properties not specific to the isotope (e.g., *x-ray scattering factors*)
    are retrieved from the associated element.
    """
    element: "Element"
    ion: IonSet # TODO: should be IonSet["Isotope"]
    isotope: int

    # mass.py
    # TODO: Do we still need to modify isotope.abundance during mass.init()?
    @property
    def abundance(self) -> float: return self._abundance
    _abundance: float
    _abundance_unc: float
    abundance_units: str # isotope only

    # activation.py
    neutron_activation: tuple["ActivationResult"]  # isotope only

    # nsf.py
    nuclear_spin: str

    def __init__(self, element: "Element", isotope_number: int):
        self.element = element
        self.isotope = isotope_number
        self.ion = IonSet(self)
    def __getattr__(self, attr):
        return getattr(self.element, attr)
    def __str__(self) -> str:
        # Deuterium and Tritium are special
        if 'symbol' in self.__dict__:
            return self.symbol
        return "%d-%s"%(self.isotope, self.element.symbol)
    def __repr__(self) -> str:
        return "%s[%d]"%(self.element.symbol, self.isotope)
    def __reduce__(self):
        return _make_isotope, (self.element.table,
                               self.element.number,
                               self.isotope)

class Element(_AtomBase):
    """
    Periodic table entry for an element.

    An element is a name, symbol and number, plus a set of properties.
    Individual isotopes can be referenced as element[*isotope_number*].
    Individual ionization states can be referenced by element.ion[*charge*].
    """
    table: str = PUBLIC_TABLE_NAME
    charge: int = 0

    def __init__(self, name: str, symbol: str, Z: int, ions: tuple[int, ...], table: "str"):
        self.name = name
        self.symbol = symbol
        self.number = Z
        self._isotopes: dict[int, "Isotope"] = {} # The actual isotopes
        self.ions = ions
        self.ion = IonSet(self)
        # Remember the table name for pickle dump/load
        if table != self.table:
            self.table = table

    @property
    def isotopes(self) -> list[int]:
        """List of all isotopes"""
        # Note: may want to return the iterator rather than the list...
        return list(sorted(self._isotopes.keys()))

    def add_isotope(self, number: int) -> "Isotope":
        """
        Add an isotope for the element.

        :Parameters:
            *number* : integer
                Isotope number, which is the number protons plus neutrons.

        :Returns: None
        """
        if number not in self._isotopes:
            self._isotopes[number] = Isotope(self, number)
        return self._isotopes[number]

    def __getitem__(self, number: int) -> "Isotope":
        try:
            return self._isotopes[number]
        except KeyError:
            raise KeyError("%s is not an isotope of %s"%(number, self.symbol))

    def __iter__(self) -> Iterator["Isotope"]:
        """
        Process the isotopes in order
        """
        for _, iso in sorted(self._isotopes.items()):
            yield iso

    # Note: using repr rather than str for the element symbol so
    # that lists of elements print nicely.  Since elements are
    # effectively singletons, the symbol name is the representation
    # of the instance.
    def __repr__(self) -> str:
        return self.symbol

    def __reduce__(self):
        return _make_element, (self.table, self.number)

def isatom(val: Any) -> bool:
    """Return true if value is an element, isotope or ion"""
    return isinstance(val, (Element, Isotope, Ion))

def isisotope(val: Any) -> bool:
    """Return true if value is an isotope or isotope ion."""
    if ision(val):
        val = val.element
    return isinstance(val, Isotope)

def ision(val: Any) -> bool:
    """Return true if value is a specific ion of an element or isotope"""
    return isinstance(val, Ion)

def iselement(val: Any) -> bool:
    """Return true if value is an element or ion in natural abundance"""
    if ision(val):
        val = val.element
    return isinstance(val, Element)

def change_table(atom: AtomVar, table: PeriodicTable) -> Atom:
    """Search for the same element, isotope or ion from a different table"""
    if ision(atom):
        if isisotope(atom):
            iso = cast(Isotope, atom)
            return table[iso.number][iso.isotope].ion[iso.charge]
        else:
            return table[atom.number].ion[atom.charge]
    else:
        if isisotope(atom):
            iso = cast(Isotope, atom)
            return table[iso.number][iso.isotope]
        else:
            return table[atom.number]

PRIVATE_TABLES: dict[str, PeriodicTable] = {}
def _get_table(name: str) -> PeriodicTable:
    try:
        return PRIVATE_TABLES[name]
    except KeyError:
        raise ValueError("Periodic table '%s' is not initialized"%name)

def _make_element(table: str, Z: int) -> Element:
    return _get_table(table)[Z]
def _make_isotope(table: str, Z: int, n: int) -> Isotope:
    return _get_table(table)[Z][n]
def _make_ion(table: str, Z: int, c: int) -> Ion:
    return _get_table(table)[Z].ion[c]
def _make_isotope_ion(table: str, Z: int, n: int, c: int) -> Ion:
    return _get_table(table)[Z][n].ion[c]


# pylint: disable=bad-whitespace
element_base: dict[int, tuple[str, str, list[int], list[int]]] = {
    # number: name symbol common_ions uncommon_ions
    # ion info comes from Wikipedia: list of oxidation states of the elements.
    0: ('neutron',     'n',  [],         []),
    1: ('Hydrogen',    'H',  [-1, 1],    []),
    2: ('Helium',      'He', [],         [1, 2]),  # +1,+2  http://periodic.lanl.gov/2.shtml
    3: ('Lithium',     'Li', [1],        []),
    4: ('Beryllium',   'Be', [2],        [1]),
    5: ('Boron',       'B',  [3],        [-5, -1, 1, 2]),
    6: ('Carbon',      'C',  [-4, -3, -2, -1, 1, 2, 3, 4], []),
    7: ('Nitrogen',    'N',  [-3, 3, 5], [-2, -1, 1, 2, 4]),
    8: ('Oxygen',      'O',  [-2],       [-1, 1, 2]),
    9: ('Fluorine',    'F',  [-1],       []),
    10: ('Neon',       'Ne', [],         []),
    11: ('Sodium',     'Na', [1],        [-1]),
    12: ('Magnesium',  'Mg', [2],        [1]),
    13: ('Aluminum',   'Al', [3],        [-2, -1, 1, 2]),
    14: ('Silicon',    'Si', [-4, 4],    [-3, -2, -1, 1, 2, 3]),
    15: ('Phosphorus', 'P',  [-3, 3, 5], [-2, -1, 1, 2, 4]),
    16: ('Sulfur',     'S',  [-2, 2, 4, 6],    [-1, 1, 3, 5]),
    17: ('Chlorine',   'Cl', [-1, 1, 3, 5, 7], [2, 4, 6]),
    18: ('Argon',      'Ar', [],         []),
    19: ('Potassium',  'K',  [1],        [-1]),
    20: ('Calcium',    'Ca', [2],        [1]),
    21: ('Scandium',   'Sc', [3],        [1, 2]),
    22: ('Titanium',   'Ti', [4],        [-2, -1, 1, 2, 3]),
    23: ('Vanadium',   'V',  [5],        [-3, -1, 1, 2, 3, 4]),
    24: ('Chromium',   'Cr', [3, 6],     [-4, -2, -1, 1, 2, 4, 5]),
    25: ('Manganese',  'Mn', [2, 4, 7],  [-3, -2, -1, 1, 3, 5, 6]),
    26: ('Iron',       'Fe', [2, 3, 6],  [-4, -2, -1, 1, 4, 5, 7]),
    27: ('Cobalt',     'Co', [2, 3],     [-3, -1, 1, 4, 5]),
    28: ('Nickel',     'Ni', [2],        [-2, -1, 1, 3, 4]),
    29: ('Copper',     'Cu', [2],        [-2, 1, 3, 4]),
    30: ('Zinc',       'Zn', [2],        [-2, 1]),
    31: ('Gallium',    'Ga', [3],        [-5, -4, -2, -1, 1, 2]),
    32: ('Germanium',  'Ge', [-4, 2, 4], [-3, -2, -1, 1, 3]),
    33: ('Arsenic',    'As', [-3, 3, 5], [-2, -1, 1, 2, 4]),
    34: ('Selenium',   'Se', [-2, 2, 4, 6], [-1, 1, 3, 5]),
    35: ('Bromine',    'Br', [-1, 1, 3, 5], [4, 7]),
    36: ('Krypton',    'Kr', [2],        []),
    37: ('Rubidium',   'Rb', [1],        [-1]),
    38: ('Strontium',  'Sr', [2],        [1]),
    39: ('Yttrium',    'Y',  [3],        [1, 2]),
    40: ('Zirconium',  'Zr', [4],        [-2, 1, 2, 3]),
    41: ('Niobium',    'Nb', [5],        [-3, -1, 1, 2, 3, 4]),
    42: ('Molybdenum', 'Mo', [4, 6],     [-4, -2, -1, 1, 2, 3, 5]),
    43: ('Technetium', 'Tc', [4, 7],     [-3, -1, 1, 2, 3, 5, 6]),
    44: ('Ruthenium',  'Ru', [3, 4],     [-4, -2, 1, 2, 5, 6, 7, 8]),
    45: ('Rhodium',    'Rh', [3],        [-3, -1, 1, 2, 4, 5, 6]),
    46: ('Palladium',  'Pd', [2, 4],     [1, 3, 5, 6]),
    47: ('Silver',     'Ag', [1],        [-2, -1, 2, 3, 4]),
    48: ('Cadmium',    'Cd', [2],        [-2, 1]),
    49: ('Indium',     'In', [3],        [-5, -2, -1, 1, 2]),
    50: ('Tin',        'Sn', [-4, 2, 4], [-3, -2, -1, 1, 3]),
    51: ('Antimony',   'Sb', [-3, 3, 5], [-2, -1, 1, 2, 4]),
    52: ('Tellurium',  'Te', [-2, 2, 4, 6], [-1, 1, 3, 5]),
    53: ('Iodine',     'I',  [-1, 1, 3, 5, 7], [4, 6]),
    54: ('Xenon',      'Xe', [2, 4, 6],  [8]),
    55: ('Cesium',     'Cs', [1],        [-1]),
    56: ('Barium',     'Ba', [2],        [1]),
    57: ('Lanthanum',  'La', [3],        [1, 2]),
    58: ('Cerium',     'Ce', [3, 4],     [2]),
    59: ('Praseodymium', 'Pr', [3],      [2, 4, 5]),
    60: ('Neodymium',  'Nd', [3],        [2, 4]),
    61: ('Promethium', 'Pm', [3],        [2]),
    62: ('Samarium',   'Sm', [3],        [2]),
    63: ('Europium',   'Eu', [2, 3],     []),
    64: ('Gadolinium', 'Gd', [3],        [1, 2]),
    65: ('Terbium',    'Tb', [3],        [1, 2, 4]),
    66: ('Dysprosium', 'Dy', [3],        [2, 4]),
    67: ('Holmium',    'Ho', [3],        [2]),
    68: ('Erbium',     'Er', [3],        [2]),
    69: ('Thulium',    'Tm', [3],        [2]),
    70: ('Ytterbium',  'Yb', [3],        [2]),
    71: ('Lutetium',   'Lu', [3],        [2]),
    72: ('Hafnium',    'Hf', [4],        [-2, 1, 2, 3]),
    73: ('Tantalum',   'Ta', [5],        [-3, -1, 1, 2, 3, 4]),
    74: ('Tungsten',   'W',  [4, 6],     [-4, -2, -1, 1, 2, 3, 5]),
    75: ('Rhenium',    'Re', [4],        [-3, -1, 1, 2, 3, 5, 6, 7]),
    76: ('Osmium',     'Os', [4],        [-4, -2, -1, 1, 2, 3, 5, 6, 7, 8]),
    77: ('Iridium',    'Ir', [3, 4],     [-3, -1, 1, 2, 5, 6, 7, 8, 9]),
    78: ('Platinum',   'Pt', [2, 4],     [-3, -2, -1, 1, 3, 5, 6]),
    79: ('Gold',       'Au', [3],        [-3, -2, -1, 1, 2, 5]),
    80: ('Mercury',    'Hg', [1, 2],     [-2, 4]),  # +4  doi:10.1002/anie.200703710
    81: ('Thallium',   'Tl', [1, 3],     [-5, -2, -1, 2]),
    82: ('Lead',       'Pb', [2, 4],     [-4, -2, -1, 1, 3]),
    83: ('Bismuth',    'Bi', [3],        [-3, -2, -1, 1, 2, 4, 5]),
    84: ('Polonium',   'Po', [-2, 2, 4], [5, 6]),
    85: ('Astatine',   'At', [-1, 1],    [3, 5, 7]),
    86: ('Radon',      'Rn', [2],        [6]),
    87: ('Francium',   'Fr', [1],        []),
    88: ('Radium',     'Ra', [2],        []),
    89: ('Actinium',   'Ac', [3],        []),
    90: ('Thorium',        'Th', [4],    [1, 2, 3]),
    91: ('Protactinium',   'Pa', [5],    [3, 4]),
    92: ('Uranium',        'U',  [6],    [1, 2, 3, 4, 5]),
    93: ('Neptunium',      'Np', [5],    [2, 3, 4, 6, 7]),
    94: ('Plutonium',      'Pu', [4],    [2, 3, 5, 6, 7]),
    95: ('Americium',      'Am', [3],    [2, 4, 5, 6, 7]),
    96: ('Curium',         'Cm', [3],    [4, 6]),
    97: ('Berkelium',      'Bk', [3],    [4]),
    98: ('Californium',    'Cf', [3],    [2, 4]),
    99: ('Einsteinium',    'Es', [3],    [2, 4]),
    100: ('Fermium',       'Fm', [3],    [2]),
    101: ('Mendelevium',   'Md', [3],    [2]),
    102: ('Nobelium',      'No', [2],    [3]),
    103: ('Lawrencium',    'Lr', [3],    []),
    104: ('Rutherfordium', 'Rf', [4],    []),
    105: ('Dubnium',       'Db', [5],    []),
    106: ('Seaborgium',    'Sg', [6],    []),
    107: ('Bohrium',       'Bh', [7],    []),
    108: ('Hassium',       'Hs', [8],    []),
    109: ('Meitnerium',    'Mt', [],     []),
    110: ('Darmstadtium',  'Ds', [],     []),
    111: ('Roentgenium',   'Rg', [],     []),
    112: ('Copernicium',   'Cn', [2],    []),
    113: ('Nihonium',      'Nh', [],     []),
    114: ('Flerovium',     'Fl', [],     []),
    115: ('Moscovium',     'Mc', [],     []),
    116: ('Livermorium',   'Lv', [],     []),
    117: ('Tennessine',    'Ts', [],     []),
    118: ('Oganesson',     'Og', [],     []),
}
# pylint: enable=bad-whitespace

def default_table(table: PeriodicTable|None=None) -> PeriodicTable:
    """
    Return the default table unless a specific table has been requested.

    This is to be used in a context like::

        def summary(table=None):
            table = core.default_table(table)
            ...
    """
    return table if table is not None else PUBLIC_TABLE

def define_elements(table: PeriodicTable, namespace: dict[str, Any]) -> list[str]:
    """
    Define external variables for each element in namespace. Elements
    are defined both by name and by symbol.

    This is called from *__init__* as::

        elements = core.default_table()
        __all__  += core.define_elements(elements, globals())

    :Parameters:
         *table* : PeriodicTable
             Set of elements
         *namespace* : dict
             Namespace in which to add the symbols.
    :Returns: [string, ...]
        A sequence listing the names defined.

    .. Note:: This will only work for *namespace* globals(), not locals()!
    """

    # Build the dictionary of element symbols
    names: list[str] = []
    for atom in [*table, table.D, table.T, table.n]:
        namespace[atom.symbol] = atom
        namespace[atom.name] = atom
        names.extend((atom.symbol, atom.name))
    return names

def get_data_path(data: Path|str) -> str:
    """
    Locate the directory for the tables for the named extension.

    :Parameters:
         *data* : string
              Name of the extension data directory.  For example, the xsf
              extension has data in the 'xsf' data directory.

    :Returns: string
         Path to the data.
    """
    import sys
    import os

    # Check for data path in the environment
    key = 'PERIODICTABLE_DATA'
    if key in os.environ:
        path = os.path.join(os.environ[key], data)
        if not os.path.isdir(path):
            raise RuntimeError('Path in environment %s not a directory'%key)
        return path

    # Check for data path in the package
    path = os.path.join(os.path.dirname(__file__), data)
    if os.path.isdir(path):
        return path

    # Check for data path next to exe/zip file.
    exepath = os.path.dirname(sys.executable)
    path = os.path.join(exepath, 'periodictable-data', data)
    if os.path.isdir(path):
        return path

    # py2app puts the data in Contents/Resources, but the executable
    # is in Contents/MacOS.
    path = os.path.join(exepath, '..', 'Resources', 'periodictable-data', data)
    if os.path.isdir(path):
        return path

    raise RuntimeError('Could not find the periodic table data files')

# Make a common copy of the table for everyone to use --- equivalent to
# a singleton without incurring any complexity.
PUBLIC_TABLE: PeriodicTable = PeriodicTable(PUBLIC_TABLE_NAME)
