import periodictable
from periodictable.core import Element, Isotope

# Test symbol access
c_element: Element = periodictable.C
fe_element: Element = periodictable.Fe
u_element: Element = periodictable.U

# Test name access
carbon: Element = periodictable.carbon
iron: Element = periodictable.iron
uranium: Element = periodictable.uranium

# Test that we can access element properties
mass: float = periodictable.C.mass
symbol: str = periodictable.Fe.symbol
number: int = periodictable.U.number

# Test special cases
deuterium: Isotope = periodictable.D
tritium: Isotope = periodictable.T
neutron: Element = periodictable.n

# Test isotopes have isotope attribute
D_iso: int = periodictable.D.isotope
Dion_iso: int = periodictable.D.ion[1].isotope
Fe56_iso: int = periodictable.Fe[56].isotope

# Test that elements do not
assert not hasattr(periodictable.H, "isotope")
assert not hasattr(periodictable.H.ion[1], "isotope")
