# This program is public domain
# Author: Paul Kienzle
r"""

The following properties are added:

*    density, density_units (|g/cm^3|)
        Densities for solids and liquids are given as specific
        gravities at 20\ |deg| C unless other wise indicated by
        *density_caveat*. Densities for gaseous elements
        are given for the liquids at their boiling points.
        Missing data are represented by *None*.
*   density_caveat
        Comments on the density, if not taken in standard conditions.
*   interatomic_distance, interatomic_distance_units (|Ang|)
        Interatomic distance estimated from element density.
*   number_density, number_density_units (|1/cm^3|)
        Number density estimated from mass and density.

Density for the isotope is computed assuming that the atomic spacing
is the same as that for the element in the natural abundance.

.. doctest::

    >>> from periodictable import D, H
    >>> print("H: %.4f, D: %.4f"%(H.density, D.density))
    H: 0.0708, D: 0.1415
    >>> print(f"{(D.density/H.density) / (D.mass/H.mass):.4f}")
    1.0000

The following plot shows density for all elements:

.. plot:: plots/density_plot.py

From the X-ray data book: http://xdb.lbl.gov/Section5/Sec_5-2.html

Data were taken mostly from [#Lide1999]_. These values are reproduced in [#ILL]_.

.. [#Lide1999] Lide. D. R., Ed., CRC Handbook of Chemistry and Physics, 80th ed.
       (CRC Press, Boca Raton, Florida, 1999)
.. [#ILL] The ILL Neutron Data Booklet, Second Edition.
"""

from typing import cast

from .core import Element, Isotope, Atom, PeriodicTable, isisotope, ision
from .constants import avogadro_number

def density(atom: Atom) -> float:

    """
    Element density for natural abundance. For isotopes, return
    the equivalent density assuming identical inter-atomic spacing as the
    naturally occuring material.

    :Parameters:
        *iso_el* : isotope or element
            Name of the element or isotope.

    :Returns:
        *density* : float | |g/cm^3|

    Reference:
        *ILL Neutron Data Booklet, original values from CRC Handbook of Chemistry and Physics,
        80th ed. (1999).*

    """
    # Note: materials with pure isotopes are adjusted by the natural mass ratio.
    if isisotope(atom):
        iso = cast(Isotope, atom)
        return iso.element._density * (iso.mass/iso.element.mass)
    return cast(Element, atom)._density

def interatomic_distance(element: Element) -> float|None:
    r"""
    Estimated interatomic distance from atomic weight and density. The
    distance between isotopes is assumed to match that between atoms in
    the natural abundance.

    :Parameters:
        *element* : Element
            The element whose interatomic distance is to be calculated.

    :Returns:
        *distance* : float | |Ang|
            Estimated interatomic distance.

    Interatomic distance is computed using:

    .. math::

        d = (m/(\rho_m N_A 10^{-24}))^{1/3}

    with units:

    .. math::

        (\rm (g\cdot mol^{-1})
        / ( (g\cdot cm^{-3})
            (atoms\cdot mol^{-1})
            (10^{-8} cm\cdot \AA^{-1})^3))^{1/3} = \AA

    """
    # Note: don't need to check for Ion or Isotope because these delegate to Element
    if element.density is None or element.mass is None:
        return None
    return (element.mass/(element.density*avogadro_number*1e-24))**(1./3.)

def number_density(element: Element) -> float:
    r"""
    Estimate the number density from atomic weight and density. The density
    for isotopes is assumed to match that of between atoms in natural abundance.

    :Parameters:
        *element* : element
            Name of the element whose number density needs to be calculated.

    :Returns:
        *Nb* : float | |1/cm^3|
            Number density of a element.

    Number density is computed using:

    .. math::

        d = N_A \frac{\rho}{m}

    with units:

    .. math::

        \rm (atoms\cdot mol^{-1})  (g\cdot cm^{-3}) / (g\cdot mol^{-1})
            = atoms\cdot cm^{-3}

    """
    # Note: don't need to check for Ion or Isotope because these delegate to Element
    if element.density is None or element.mass is None:
        return None
    return (element.density/element.mass)*avogadro_number

def init(table: PeriodicTable, reload: bool=False) -> None:
    if 'density' in table.properties and not reload:
        return
    table.properties.append('density')
    Isotope.density \
        = property(density, doc="density using inter-atomic spacing from naturally occurring form")  # type: ignore[assignment]
    Element.density \
        = property(density, doc="density using inter-atomic spacing from naturally occurring form")  # type: ignore[assignment]
    Element.density_units = "g/cm^3"

    Element.interatomic_distance \
        = property(interatomic_distance,
                   doc="interatomic distance estimated from density")  # type: ignore[assignment]
    Element.interatomic_distance_units = "angstrom"
    Element.number_density \
        = property(number_density,
                   doc="number density estimated from mass and density")  # type: ignore[assignment]
    Element.number_density_units = "1/cm^3"

    for k, v in element_densities.items():
        el = getattr(table, k)
        if isinstance(v, tuple):
            el._density = v[0]
            el.density_caveat = v[1]
        elif v is None:
            el._density = None
            el.density_caveat = "unavailable"
        else:
            el._density = v
            el.density_caveat = ""

element_densities = dict(
    n=None,  # No longer using element 0 for a bare neutron
    H=(0.0708, "T=-252.87"),
    He=(0.122, "T=-268.93"),
    Li=0.534,
    Be=1.848,
    B=2.34,
    C=(2.2, "1.9-2.3 (graphite); CXRO, RSC and CRC use 2.2; NIST Xray tables use 2.26"),
    N=(0.808, "T=-195.79"),
    O=(1.14, "T=-182.95"),
    F=(1.50, "T=-188.12"),
    Ne=(1.207, "T=-246.08"),
    Na=0.971,
    Mg=1.738,
    Al=2.6989,
    Si=(2.33, "T=25"),
    P=1.82,
    S=2.07,
    Cl=(1.56, "T=-33.6, 0.44 C above boiling point"),
    Ar=(1.40, "T=-185.85"),
    K=0.862,
    Ca=1.55,
    Sc=(2.989, "T=25"),
    Ti=4.54,
    V=(6.11, "T=18.7"),
    Cr=(7.19, "7.18-7.20"),
    Mn=(7.33, "7.21-7.44"),
    Fe=7.874,
    Co=8.9,
    Ni=(8.902, "T=25"),
    Cu=8.96,
    Zn=(7.133, "T=25"),
    Ga=(5.904, "T=29.6"),
    Ge=(5.323, "T=25"),
    As=5.73,
    Se=4.79,
    Br=3.12,
    Kr=(2.16, "T=-153.22"),
    Rb=1.532,
    Sr=2.54,
    Y=(4.469, "T=25"),
    Zr=6.506,
    Nb=8.57,
    Mo=10.22,
    Tc=(11.50, "calculated"),
    Ru=12.41,
    Rh=12.41,
    Pd=12.02,
    Ag=10.50,
    Cd=8.65,
    In=7.31,
    Sn=7.31,
    Sb=6.691,
    Te=6.24,
    I=4.93,
    Xe=(3.52, "T=-108.12"),
    Cs=1.873,
    Ba=3.5,
    La=(6.145, "T=25"),
    Ce=(6.770, "T=25"),
    Pr=6.773,
    Nd=(7.008, "T=25"),
    Pm=(7.264, "T=25"),
    Sm=(7.520, "T=25"),
    Eu=(5.244, "T=25"),
    Gd=(7.901, "T=25"),
    Tb=8.230,
    Dy=(8.551, "T=25"),
    Ho=(8.795, "T=25"),
    Er=(9.066, "T=25"),
    Tm=(9.321, "T=25"),
    Yb=6.966,
    Lu=(9.841, "T=25"),
    Hf=13.31,
    Ta=16.654,
    W=19.3,
    Re=21.02,
    Os=22.57,
    Ir=(22.42, "T=17"),
    Pt=21.45,
    Au=(19.3, "approximate"),
    Hg=13.546,
    Tl=11.85,
    Pb=11.35,
    Bi=9.747,
    Po=9.32,
    At=None,
    Rn=None,
    Fr=None,
    Ra=None,
    Ac=None,
    Th=11.72,
    Pa=(15.37, "calculated"),
    U=(18.95, "approximate"),
    Np=20.25,
    Pu=(19.84, "T=25"),
    Am=13.67,
    Cm=(13.51, "calculated"),
    Bk=(14, "estimated"),
    Cf=None,
    Es=None,
    Fm=None,
    Md=None,
    No=None,
    Lr=None,
    Rf=None,
    Db=None,
    Sg=None,
    Bh=None,
    Hs=None,
    Mt=None,
    Ds=None,
    Rg=None,
    Cn=None,
    Nh=None,
    Fl=None,
    Mc=None,
    Lv=None,
    Ts=None,
    Og=None,
    )
