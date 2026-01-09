# -*- coding: utf-8 -*-
# This program is public domain
# Author: Paul Kienzle
r"""

Neutron scattering factors for the elements and isotopes.

For details of neutron scattering factor values, see :class:`Neutron`.
The property is set to *None* if there is no neutron scattering information
for the element. Individual isotopes may have their own scattering
information.

Example
=======

Print a table of coherent scattering length densities for isotopes
of a particular element:

.. doctest::

    >>> import periodictable
    >>> for iso in periodictable.Ni:
    ...     if iso.neutron.has_sld():
    ...         print("%s %7.4f"%(iso,iso.neutron.sld()[0]))
    58-Ni 13.1526
    60-Ni  2.5575
    61-Ni  6.9417
    62-Ni -7.9464
    64-Ni -0.3379


Details
=======

There are a number of functions available in periodictable.nsf

    :func:`neutron_energy`
        Return neutron energy given wavelength.

    :func:`neutron_wavelength`
        Return wavelength given neutron energy.

    :func:`neutron_wavelength_from_velocity`
        Return wavelength given neutron velocity.

    :func:`neutron_scattering`
        Computes scattering length density, cross sections and
        penetration depth for a compound.

    :func:`neutron_sld`
        Computes scattering length density for a compound.

    :func:`neutron_composite_sld`
        Returns a scattering length density for a compound whose composition
        is variable.

    :func:`energy_dependent_table`
        Lists isotopes with energy dependence.

    :func:`sld_table`
        Lists scattering length densitys for all elements in natural abundance.

    :func:`absorption_comparison_table`
        Compares the imaginary bound coherent scattering length to the
        absorption cross section.

    :func:`coherent_comparison_table`
        Compares the bound coherent scattering length to the
        coherent scattering cross section.

    :func:`total_comparison_table`
        Compares the total scattering cross section to the sum of the
        coherent and incoherent scattering cross sections.

For private tables use :func:`init` to set the data.

The neutron scattering information table is reproduced from the Atomic
Institute for Austrian Universities\ [#Rauch2003]_,\ [#Rauch2000]_:

    `<http://www.ati.ac.at/~neutropt/scattering/table.html>`_

The above site has references to the published values for every entry in
the table.  We have included these in the documentation directory
within the periodictable source package. Some typographical errors have
been fixed. In particular, Zn-70 has b_c listed as 6.9 in the table,
but 6.0 in the source materials for the table.

Note that enteries in the table have been measured independently, so the values
for the scattering length of an element or isotope may be inconsistent with
the values measured for the corresponding cross section. The comparison table
functions highlight these differences.

Tables from Sears\ [#Sears1992]_\ [#Sears2006]_, Rauch\ [#Rauch2003]_ and
Dawidowski\ [#Dawidowski2013]_ make different choices for the recommended
values. These are noted in periodictable issue #59
`<https://github.com/python-periodictable/periodictable/issues/59>`_
with changes from Sears to Rauch
`(a) <https://github.com/python-periodictable/periodictable/issues/59#issuecomment-1693686953>`__
and from Rauch to Dawidowski
`(b) <https://github.com/python-periodictable/periodictable/issues/59#issuecomment-1690212205>`__.

The following newer measurements from the literature are included:

    1H b_c -3.7423(12) => -3.7395(11) [1]
    2H b_c 6.674(6) => 6.6681(27) [1]
    4He b_c 3.26(3) => 3.0982(21) [3] (see also [5], which gives 3.075(6))
    4He coherent = total cross sections computed from 4 pi b_c^2/100
    natHe computed from isotopic weighting of 3He and 4He
    natC b_c 6.6484(13) => 6.6472(9) [1]
    natO b_c 5.805(4) => 5.8037(29) [1]
    17O b_c 5.6(5) => 5.867(4) [2]
    18O b_c 5.84(7) => 6.009(5) [2]
    natSn b_c 6.225(2) => 6.2239(13) [1]
    154Sm b_c 8.0(1.0) => 8.97(6) [4]
    153Eu b_c 8.22(12) => 8.85(3) [4]
    191Ir b_c => 12.1(9) [6]
    193Ir b_c => 9.71(18) [6]
    natPb b_c 9.401(2) => 9.4024(13) [1]
    natBi b_c 8.532(2) => 8.5242(18) [1]

    He total 1.34(2) => 1.188(5) [5] (ignored; using 4 pi b_c^2/100)
    Ar total 0.683(4) => 0.683(5) [5] (ignored; existing value is more precise)
    Kr total 7.68(13) => 7.685(26) [5]
    Xe total --- => 4.344(17) [5]

    [1] Snow (2020) 10.1103/PhysRevD.101.062004 [average of values in Table 1]
    [2] Fischer (2012) 10.1088/0953-8984/24/50/505105
    [3] Haun (2020) 10.1103/PhysRevLett.124.012501
    [4] Kohlmann (2016) 10.1515/zkri-2016-1984
    [5] Haddock (2019) 10.1103/PhysRevC.100.064002
    [6] Hannon (2018) 10.1107/S1600576718006064

.. [#Rauch2003] Rauch, H. and Waschkowski, W. (2003)
    Neutron Scattering Lengths in ILL
    Neutron Data Booklet (second edition), A.-J. Dianoux, G. Lander, Eds.
    Old City Publishing, Philidelphia, PA. pp 1.1-1 to 1.1-17.
    (https://www.ill.eu/fileadmin/user_upload/ILL/1_About_ILL/Documentation/NeutronDataBooklet.pdf
    Retrieved March 2008)

.. [#Rauch2000] Rauch, H. and Waschkowski, W. (2000)
    Neutron scattering lengths. Schopper, H. (ed.). SpringerMaterials -
    The Landolt-Börnstein Database (http://www.springermaterials.com).
    doi:10.1007/10499706_6

.. [#Lynn1990] Lynn, J.E. and Seeger, P.A. (1990)
    Resonance effects in neutron scattering lengths of rare-earth nuclides.
    Atomic Data and Nuclear Data Tables 44, 191-207.
    doi:10.1016/0092-640X(90)90013-A

.. [#Sears2006] Sears, V. F. (2006)
    4.4.4 Scattering lengths for neutrons.
    In Prince, E. ed. Intl. Tables for Crystallography C.
    Kluwer Academic Publishers. pp 444-454.
    (https://it.iucr.org/Cb/ch4o4v0001/sec4o4o4/)
    doi:10.1107/97809553602060000103

.. [#Sears1992] Sears, V.F. (1992)
    Neutron scattering lengths and cross sections.
    Neutron News 3, No. 3, 26-37.

.. [#May1982] May, R.P.,  Ibel, K. and Haas, J. (1982)
    The forward scattering of cold neutrons by mixtures of light and heavy water.
    J. Appl. Cryst. 15, 15-19. doi:10.1107/S0021889882011285

.. [#Mildner1998] Mildner, D.F.R., Lamaze, G.P. (1998)
   Neutron Transmission of Single-Crystal Sapphire.
   J Appl Crystallogr 31, 835–840. doi:10.1107/S0021889898005846

.. [#Glinka2011] Glinka, C.J. (2011)
    Incoherent Neutron Scattering from Multi-element Materials.
    J. Appl. Cryst. 44, 618-624. doi:10.1107/S0021889811008223

.. [#Dawidowski2013] Dawidowski, J., Granada, J. R., Santisteban,
    J. R., Cantargi, F., & Palomino, L. A. R. (2013).
    Appendix—Neutron Scattering Lengths and Cross Sections.
    In F. Fernandez-Alonso & D. L. Price (Eds.),
    Experimental Methods in the Physical Sciences (Vol. 44, pp. 471–528).
    Academic Press. doi:10.1016/B978-0-12-398374-9.09989-7
"""

#.. [#Koester1991] Koester, L., Rauch, H., Seymann. E. (1991)
#    Atomic Data Nuclear Data Tables 49, 65. doi:10.1016/0092-640X(91)90012-S
#.. [#Smith2006] Smith, G.S. and Majkrzak, C.M. (2006)
#    2.9 Neutron reflectometry.
#    In E. Prince ed. Intl. Tables for Crystallography C.
#    Wiley InterScience. pp 126-146. doi:10.1107/97809553602060000584
#

from typing import TYPE_CHECKING, cast

import numpy as np
from numpy import sqrt, pi, asarray, inf
from numpy.typing import NDArray, ArrayLike
from .core import Element, Isotope, PeriodicTable, default_table
from .constants import (avogadro_number, planck_constant, electron_volt,
                        neutron_mass, atomic_mass_constant)

DoubleArray = NDArray[np.float64]

__all__ = ['init', 'Neutron',
           'neutron_energy', 'neutron_wavelength',
           'neutron_wavelength_from_velocity',
           'neutron_scattering', 'neutron_sld', 'neutron_composite_sld',
           'sld_plot',
           'absorption_comparison_table', 'coherent_comparison_table',
           'incoherent_comparison_table', 'total_comparison_table',
           'energy_dependent_table', 'sld_table',
           'neutron_sld_from_atoms',
           #'scattering_potential',
          ]

#: Wavelength [Å] for which neutron scattering cross sections are tabulated.
ABSORPTION_WAVELENGTH = 1.798 # [Å]


#: Energy [eV] <=> wavelength [Å]:
#:   E = 1/2 m v² = h² / (2 m λ²)
#:   E[meV s]  = h[J s]²/(2 m_n[u] m_u[kg/u] λ[Å]²/10^20[Å/m])
#:             * 1000[meV/eV] / electron_volt[J/eV]
ENERGY_FACTOR = (
    1e23 * planck_constant**2/electron_volt
    / (2 * neutron_mass * atomic_mass_constant))

#: Velocity[m/s] <=> wavelength[Å]:
#:  h = p λ = m v λ
#:  λ[Å] = h[J s] / ( m_n[kg] v[m/s] ) 10^10[Å/m]
VELOCITY_FACTOR = (
    1e10 * planck_constant / (neutron_mass * atomic_mass_constant))

def neutron_wavelength(energy: ArrayLike) -> NDArray:
    r"""
    Convert neutron energy to wavelength.

    :Parameters:
        *energy* \: float or vector | meV

    :Returns:
        *wavelength* \: float or vector | |Ang|

    Energy is converted to wavelength using

    .. math::

        E = 1/2 m_n v^2 = h^2 / (2 m_n \lambda^2)
        \Rightarrow \lambda = \sqrt{h^2 / (2 m_n E)}

    where

        $h$ = Planck constant in |Js|

        $m_n$ = neutron mass in kg

    """
    return sqrt(ENERGY_FACTOR / asarray(energy))

# TODO: why is neutron_wavelength_from_velocity not using asarray() ?
def neutron_wavelength_from_velocity(velocity: float) -> float:
    r"""
    Convert neutron velocity to wavelength.

    :Parameters:
        *velocity* \: float or vector | m/s

    :Returns:
        *wavelength* \: float or vector | |Ang|

    Velocity is converted to wavelength using

    .. math::

        \lambda = h/p = h/(m_n v)

    where

        $h$ = Planck constant in |Js|

        $m_n$ = neutron mass in kg
    """
    return VELOCITY_FACTOR / velocity

def neutron_energy(wavelength: ArrayLike) -> NDArray:
    r"""
    Convert neutron wavelength to energy.

    :Parameters:
        *wavelength* \: float or vector | |Ang|

    :Returns:
        *energy* \: float or vector | meV

    Wavelength is converted to energy using

    .. math::

        E = 1/2 m_n v^2 = h^2 / (2 m_n \lambda^2)

    where:

        $h$ = Planck constant in |Js|

        $m_n$ = neutron mass in kg
    """
    return ENERGY_FACTOR / asarray(wavelength)**2

def _CHECK_scattering_potential(sld):
    r"""
    Convert neutron scattering length density to energy potential.

    :Parameters:
        *sld* \: float or vector | |1e-6/Ang^2|

            Scattering length density.

    :Returns:
        *energy* \: float or vector | $10^{-6}$ eV

            Scattering potential.

    Computes:[#Smith2006]_

    .. math::

        V = 2 \pi \hbar^2 N_b / m_n

    where:

        $\hbar = h / (2 \pi)$

        $h$ = Planck constant in |Js|

        $N_b = \sum{ n_i b_i } / V$

        $m_n$ = neutron mass in kg
    """
    return (ENERGY_FACTOR/pi) * asarray(sld)

_4PI_100 = 4*np.pi/100
class Neutron:
    r"""
    Neutron scattering factors are attached to each element in the periodic
    table for which values are available.  If no information is available,
    then the neutron field of the element will be *None*. Even when neutron
    information is available, it may not be complete, so individual fields
    may be *None*.

    The following fields are defined:

    * b_c (fm)
        Bounds coherent scattering length.

    * total (barn)
        Total scattering cross section $\sigma_s$.  This does not include the
        absorption cross section.  To compute the total collision cross
        section use $\sigma_t = \sigma_s + \sigma_a$

    * absorption (barn)
        Absorption cross section $\sigma_a$ at 1.798 |Ang|.  Scale to your beam
        by dividing by periodictable.nsf.ABSORPTION_WAVELENGTH and multiplying
        by your wavelength. This wavelength corresponds to a neutron velocity
        of 2200 m/s and neutron energy of 25.3 meV.

    * b_c_complex (fm)
        Complex coherent scattering length derived from the tabulated
        values using $b_c - i \sigma_a / (1000 \cdot 2 \lambda)$.

    Additional columns not used for calculation include:

    * b_c_i (fm)
        Imaginary bound coherent scattering length.  This is related to
        absorption cross section by $\sigma_a = 4 \pi \mathrm{Im}(b_c)/k$ where
        $k = 2 \pi/\lambda$ and an additional factor of 1000 for converting
        between |Ang|\ |cdot|\ fm and barns.  b_c_i is not available for
        all isotopes for which absorption cross sections have been measured.

    * bp, bm (fm)
        Spin-dependent scattering for I+1/2 and I-1/2 (not always available).
        Incoherent scattering arises from the spin-dependent scattering b+
        and b-. The Neutron Data Booklet\ [#Rauch2003]_ gives formulas for
        calculating coherent and incoherent scattering from b+ and b- alone.

    * bp_i, bm_i (fm)
        Imaginary portion of bp and bm.

    * is_energy_dependent (boolean)
        Do not use this data if scattering is energy dependent.

    * coherent (barn)
        Coherent scattering cross section.  This is tabulated but not used.
        In theory coherent scattering is related to bound coherent scattering
        by $\sigma_c = 4 \pi |\mathrm{Re}(b_c) + i \mathrm{Im}(b_c)|^2/100$.
        In practice, these values are different, with the following table
        showing the largest relative difference:

        ========  ========  ========  ========  ========
        Sc   3%   Ti   4%   V   34%   Mn   1%   Cd  2%
        Te   4%   Xe   9%   Sm  19%   Eu  44%   Tb  1%
        Ho  11%   W    4%   Au   7%   Hg   2%   Ra  3%
        ========  ========  ========  ========  ========

    * incoherent (barn)
        Incoherent scattering cross section $\sigma_i$.  This is tabulated but
        not used. Instead, the incoherent cross section is computed from the
        total cross section minus the coherent cross section even for single
        atoms so that results from compounds are consistent with results from
        single atoms.

    For elements, the scattering cross-sections are based on the natural
    abundance of the individual isotopes. Individual isotopes may have
    the following additional fields

    * abundance (%)
        Isotope abundance used to compute the properties of the element in
        natural abundance.

    * nuclear_spin (string)
        Spin on the nucleus: '0', '1/2', '3/2', etc.

    Each field ``T`` above has a corresponding ``T_units`` attribute with
    the name of the units.

    For scattering calculations the scattering length density is the value
    of interest. This is computed from the *number_density* of the individual
    elements, as derived from the element density and atomic mass.

    .. Note:: 1 barn = 100 |fm^2|
    """
    b_c: float|None = None
    b_c_units: str = "fm"
    b_c_i: float|None = None
    b_c_i_units: str = "fm"
    b_c_complex: complex|None = None
    b_c_complex_units: str = "fm"
    bp: float|None = None
    bp_i: float|None = None
    bp_units: str = "fm"
    bm: float|None = None
    bm_i: float|None = None
    bm_units: str = "fm"
    coherent: float|None = None
    coherent_units: str = "barn"
    incoherent: float|None = None
    incoherent_units: str = "barn"
    total: float|None = None
    total_units: str = "barn"
    absorption: float|None = None
    absorption_units: str = "barn"
    abundance: float = 0.
    abundance_units: str = "%"
    is_energy_dependent: bool = False
    nsf_table: tuple[DoubleArray, DoubleArray]|None = None

    def __init__(self):
        self._number_density = None
    def __str__(self) -> str:
        if not self.has_sld():
            return "unknown"
        b_c, coh, inc, abs = self.b_c, self.coherent, self.incoherent, self.absorption
        return f"b_c={cast(float, b_c):.3g} coh={cast(float, coh):.3g} inc={cast(float, inc):.3g} abs={cast(float, abs):.3g}"

    def has_sld(self) -> bool:
        """Returns *True* if sld is defined for this element/isotope."""
        # TODO: use NaN for missing information
        #return np.isnan(self.b_c * self._number_density)
        return self.b_c is not None and self._number_density is not None

    # PAK 2021-04-05: allow energy dependent b_c
    def scattering_by_wavelength(self, wavelength: ArrayLike) -> tuple[ArrayLike, ArrayLike]|tuple[None, None]:
        r"""
        Return scattering length and total cross section for each wavelength.

        For rare earth isotopes this returns the energy-dependent
        $\mathrm{Re}(b_c)$ and $\mathrm{Im}(b_c)$ interpolated into the
        scattering length tables. Values are extrapolated with constant
        values at the ends of the table. Total scattering is returned as
        $4\pi/100 |b_c|^2$ with no contribution for bound incoherent
        scattering.

        :Parameters:
            *wavelength* \: float(s) | |Ang|

        :Returns:
            *b_c* \: complex(s) | fm

            *sigma_s* \: float(s) | barn
        """
        if not self.has_sld():
            return None, None

        # TODO: do vector conversion at the end rather than the beginning.
        if self.nsf_table is None:
            ones = 1 if np.isscalar(wavelength) else np.ones_like(wavelength)
            return ones*cast(float, self.b_c_complex), ones*cast(float, self.total)
        #energy = neutron_energy(wavelength)
        #b_c = np.interp(energy, self.nsf_table[0], self.nsf_table[1])
        b_c = np.interp(wavelength, self.nsf_table[0], self.nsf_table[1])
        # TODO: sigma_s should include an incoherent contribution
        sigma_s = _4PI_100*abs(b_c)**2 # 1 barn = 1 fm^2 1e-2 barn/fm^2
        return b_c, sigma_s

    def sld(self, *, wavelength: ArrayLike=ABSORPTION_WAVELENGTH) -> NDArray|None:
        r"""
        Returns scattering length density for the element at natural
        abundance and density.

        :Parameters:
            *wavelength* \: float(s) | |Ang|

        :Returns:
            *sld* \: float(s), float(s), float(s) | |1e-6/Ang^2|
                (*real*, -*imaginary*, *incoherent*) scattering length density.

        Returns (None, None, None) if sld is not known for this element.

        See :func:`neutron_scattering` for details.
        """
        # TODO: deprecate in favour of neutron_scattering(el)
        if not self.has_sld():
            return None
        return self.scattering(wavelength=wavelength)[0]

    def scattering(self, *, wavelength: ArrayLike=ABSORPTION_WAVELENGTH) -> tuple[NDArray, NDArray, NDArray]|tuple[None, None, None]:
        r"""
        Returns neutron scattering information for the element at natural
        abundance and density.

        :Parameters:
            *wavelength* \: float(s) | |Ang|

        :Returns:
            *sld* \: float(s), float(s), float(s) | |1e-6/Ang^2|
                (*real*, -*imaginary*, *incoherent*) scattering length density

            *xs* \: float(s), float(s), float(s) | |1/cm|
                (*coherent*, *absorption*, *incoherent*) cross sections.
                :w

            *penetration* \: float(s) | cm
                1/e penetration length.

        Returns (None, None, None) if sld is not known for this element.

        See :func:`neutron_scattering` for details.
        """
        # TODO: deprecate in favour of neutron_scattering(el)
        # Compute number and absorption density assuming isotope has
        # same structure as the bulk element
        if not self.has_sld():
            return None, None, None

        number_density = self._number_density*1e-24  # N/A^3 = N/cm^3 (1e-8 cm/A)^3
        b_c, sigma_s = self.scattering_by_wavelength(wavelength)
        return _calculate_scattering(number_density, wavelength, b_c, sigma_s)

def energy_dependent_init(table: PeriodicTable) -> None:
    from .nsf_tables import ENERGY_DEPENDENT_TABLES

    for (el_name, iso_num), values in ENERGY_DEPENDENT_TABLES.items():
        energy, re_a, im_a, _ = zip(*values) # Ignoring abs(a)
        wavelength = neutron_wavelength(asarray(energy)*1000) # 1 eV = 1000 meV
        el = getattr(table, el_name)
        atom = el if iso_num is None else el[iso_num]
        xs = asarray(re_a) + 1j*asarray(im_a)
        #atom.neutron.nsf_table = asarray(energy)*1000, xs
        atom.neutron.nsf_table = wavelength[::-1], xs[::-1]
        #print(f"adding {atom}")

    # Lu nat missing from Lynn and Seeger, so mix Lu[175] and Lu[176]
    Lu175 = table.Lu[175]
    Lu176 = table.Lu[176]
    bc_175 = Lu175.neutron.b_c_complex
    wavelength, bc_176 = Lu176.neutron.nsf_table
    bc_nat = (bc_175*Lu175.abundance + bc_176*Lu176.abundance)/100.0 # 1 fm = 1fm * %/100
    table.Lu.neutron.nsf_table = wavelength, cast(DoubleArray, bc_nat)
    #table.Lu.neutron.total = 0.  # zap total cross section

def init(table: PeriodicTable, reload: bool=False) -> None:
    """
    Loads the Rauch table from the neutron data book.
    """
    if 'neutron' in table.properties and not reload:
        return
    table.properties.append('neutron')
    assert ('density' in table.properties and 'mass' in table.properties), \
        "Neutron table requires mass and density properties"

    # Defaults for missing neutron information
    missing = Neutron()
    Isotope.neutron = missing
    Element.neutron = missing

    for line in nsftable.split('\n'):
        columns = line.split(',')

        nsf = Neutron()
        p = columns[1]
        spin = columns[2]
        nsf.b_c, nsf.bp, nsf.bm = [fix_number(a) for a in columns[3:6]]
        nsf.is_energy_dependent = (columns[6] == 'E')
        nsf.coherent, nsf.incoherent, nsf.total, nsf.absorption \
            = [fix_number(a) for a in columns[7:]]
        # 1 fm = (1 barn)(100 fm^2/barn)/(1 A) (1e-5 A/fm)
        # Note: Sears (1992) uses b = b' - i b'', so negate sigma_a for b''.
        # Warning: -b_c.imag may be -0, which can mess with your calculations.
        #if nsf.b_c is None: print(f"b_c unavailable for {columns[0]}")
        b_c = nsf.b_c if nsf.b_c is not None else np.nan
        b_c_i = -nsf.absorption/(2000*ABSORPTION_WAVELENGTH)
        nsf.b_c_complex = b_c + 1j*b_c_i

        if not np.isnan(b_c):
            # Ir-191 and Ir-193 don't list scattering cross sections so deduce
            # them from the bound coherent cross section. Since these are both
            # odd isotopes there ought to be some incoherent scattering, but
            # zero is well within the uncertainty measured in bulk Ir.
            if nsf.coherent is None:
                nsf.coherent = 4*pi/100*abs(nsf.b_c_complex)**2
            if nsf.incoherent is None:
                nsf.incoherent = 0
            if nsf.total is None:
                nsf.total = nsf.coherent + nsf.incoherent

        parts = columns[0].split('-')
        Z = int(parts[0])
        symbol = parts[1]
        isotope_number = int(parts[2]) if len(parts) == 3 else 0

        # Fetch element from the table and check that the symbol matches
        element = table[Z]
        assert element.symbol == symbol, \
            "Symbol %s does not match %s" % (symbol, element.symbol)

        # Plug the default number density for the element into the nsf so
        # it can calculate sld.
        nsf._number_density = element.number_density # N/cm^3 = N/cm^3

        if isotope_number == 0:
            # Value for element using laboratory abundances of isotopes
            element.neutron = nsf
        else:
            # Values for the individual isotope
            isotope = element.add_isotope(isotope_number)
            isotope.neutron = nsf
            isotope.nuclear_spin = spin
            # p column contains either abundance(uncertainty) or "half-life Y"
            isotope.neutron.abundance = fix_number(p) if ' ' not in p else 0

            # If the element is not yet initialized, copy info into the atom.
            # This serves to set the element info for elements with only
            # one isotope.
            if element.neutron is missing:
                element.neutron = nsf

    for line in nsftableI.split('\n'):
        columns = line.split(',')

        # Fetch the nsf record
        parts = columns[0].split('-')
        Z = int(parts[0])
        symbol = parts[1]
        isotope_number = int(parts[2]) if len(parts) == 3 else 0
        element = table[Z]
        if isotope_number == 0:
            nsf = element.neutron
        else:
            nsf = element[isotope_number].neutron

        # Read imaginary values
        nsf.b_c_i, nsf.bp_i, nsf.bm_i = [
            fix_number(a) for a in columns[1:]]

    # Add energy-dependent tables
    energy_dependent_init(table)


# TODO: split incoherent into spin and isotope incoherence (eq 17-19 of Sears)
# TODO: require parsed compound rather than including formula() keywords in api
# Note: docs and function prototype are reproduced in __init__
# CRUFT: deprecated circular import with periodictable.formulas
if TYPE_CHECKING:
    from .formulas import FormulaInput
def neutron_scattering(
        compound: "FormulaInput", *,
        density: float|None=None,
        wavelength: ArrayLike|None=None,
        energy: ArrayLike|None=None,
        natural_density: float|None=None,
        table: PeriodicTable|None=None,
        ) -> tuple[NDArray, NDArray, NDArray]|tuple[None, None, None]:
    r"""
    Computes neutron scattering cross sections for molecules.

    :Parameters:
        *compound* \: Formula initializer
            Chemical formula
        *density* \: float | |g/cm^3|
            Mass density
        *natural_density* \: float | |g/cm^3|
            Mass density of formula with naturally occuring abundances
        *wavelength* 1.798 \: float(s) | |Ang|
            Neutron wavelength (default=1.798 |Ang|).
        *energy* \: float(s) | meV
            Neutron energy.  If energy is specified then wavelength is ignored.
        *table* \: PeriodicTable
            Alternate table to use when parsing *compound*.

    :Returns:
        *sld* \: float(s), float(s), float(s) | |1e-6/Ang^2|
            (*real*, -*imaginary*, *incoherent*) scattering length density.
        *xs* \: float(s), float(s), float(s) | |1/cm|
            (*coherent*, *absorption*, *incoherent*) cross sections.
        *penetration* \: float(s) | cm
            1/e penetration depth of the beam

    Returns (None, None, None) if sld is unknown for any component.

    :Raises:
        *AssertionError* \: density is missing.


    .. Note:

        The returned values will be vectors if *wavelength* is a vector.

    Neutron scattering cross sections for materials are calculated from
    tabulated values for the different nuclei. The result is only an
    approximation. Actual scattering depends on details of sample composition,
    as well as the incoming neutron energy and sample temperature, especially
    for light elements. For low energy neutrons the tabulated
    cross sections are generally a lower limit. The measured incoherent
    scattering from hydrogen, for example, can be considerably larger
    (by more than 20%) than its bound value, leading to an estimate of 5.621/cm
    for H2O as computed compared to ~7.0/cm as measured with 5 meV neutrons
    at 290K.\ [#May1982]_ The alignment of the neutron spin with the nuclei
    spin also matters, as demonstrated by $^3\mathrm{He}$ polarizers.

    The tables themselves are not self-consistent.  Because the different
    quantities are measured indirectly with a variety of techniques, there
    are discrepencies when converting values from one column to another. These
    differences can be seen with the following:

        :func:`absorption_comparison_table`

        :func:`coherent_comparison_table`

        :func:`total_comparison_table`

    For our calculations we use the real part of the bound coherent scattering
    length $\mathrm{Re}(b_c)$ (labelled b_c in the tables) and the absorption
    cross section $\sigma_a$ from which we derive the imaginary scattering
    length, $\mathrm{Im}(b_c)$. See Sears (1992) for details.\ [#Sears1992]_

    We first need to average quantities for the unit cell of the molecule.
    Molar mass *m* (g/mol) is the sum of the masses of each component:

    .. math::

        m = \sum{n_k m_k}\ {\rm for\ each\ atom}\ k=1,2,\ldots

    Cell volume $V$ (|Ang^3|/molecule) is molar mass $m$ over density
    $\rho$, with a correction based on Avogadro's number $N_A$ (atoms/mol)
    and the length conversion $10^8$ |Ang|/cm:

    .. math::

        V = m/\rho \cdot 1/N_A \cdot (10^8)^3

    Number density $N$ is the number of scatterers per unit volume:

    .. math::

        N = \left.\sum{n_k} \right/ V

    The coherent scattering length of the molecule is computed from the
    average scattering length $b_c = \mathrm{Re}(b_c) + i \mathrm{Im}(b_c)$
    weighted by frequency:

    .. math::

        \mathrm{Re}(b_c) &= \left.\sum n_k \mathrm{Re}(b_{ck}) \right/ \sum n_k \\
        \mathrm{Im}(b_c) &= \left.\sum n_k \mathrm{Im}(b_{ck}) \right/ \sum n_k

    The individual $\mathrm{Im}(b_{ck})$ values are derived from the absorption
    cross sections $\sigma_a$, tabulated at wavelength $\lambda = 1.798$ |Ang|
    and scaled to fm (with 1000 fm = 1 barn/|Ang|):

    .. math::

        \mathrm{Im}(b_{ck}) = -\left. \sigma_{ak} \right/ (1000 \cdot 2 \lambda)

    Note the sign change relative to $b''$ in Sears (1992), with Eq 2
    defining $b = b' - i b''$. Since we are not considering polarized
    nuclei, the imaginary incoherent contribution is zero and
    $b'' = -\mathrm{Im}(b_c)$.

    Some rare earth isotopes are energy-dependent with complex bound coherent
    scattering length $b_c$ tabulated by energy.\ [#Lynn1990]_ For the given
    input wavelength $\lambda$, $b_c$ is interpolated from the table values,
    with the end points used for values outside the tabulated range.

    The average scattering length is converted to scattering cross
    sections, with $\sigma_c$ scaled to barn (with 1 barn= 100 |fm^2|)
    and $\sigma_a$ scaled to barn (with 1000 barn = 1 fm |Ang|):

    .. math::

        \sigma_c &= \left. 4 \pi |\mathrm{Re}(b_c) + i \mathrm{Im}(b_c)|^2 \right/ 100 \\
        \sigma_a &= -\left. 1000 \cdot 4 \pi \left< \mathrm{Im}(b_c) \right> \right/k
            \ {\rm for} \ k=2\pi / \lambda

    For most elements the scattering length is independent of energy in
    the thermal neutron energy range so the coherent scattering length
    $\sigma_c$ is unchanged. The absorption cross section $\sigma_a$ for
    these elements scales linearly with wavelength and can be adjusted
    with a simple multiplication:

    .. math::

        \sigma_a' = \sigma_a \lambda' / \lambda_o = \sigma_a \lambda' / 1.798

    The incoherent scattering length is more complicated, including
    contributions from spin incoherence for different possible spin states
    as well as isotope incoherence from diffuse coherent
    scattering.\ [#Glinka2011]_ Using the total cross section $\sigma_s$ from
    the table, the incoherent scattering length is estimated as:

    .. math::

        \sigma_s &= \left.\sum n_k \sigma_{sk} \right/ \sum n_k \\
        \sigma_i &= \sigma_s - \sigma_c \\
        b_i &= \sqrt{ 100 \sigma_i / (4 \pi) }

    For the energy dependent rare earth isotopes the total scattering is
    estimated from $b = \mathrm{Re}(b_c) + i \mathrm{Im}(b_c)$, ignoring
    any spin incoherence effects. As a result, incoherent scattering for
    materials with energy-dependent rare earth isotopes with non-zero
    nuclear spin will be underestimated.

    The scattering potential can be expressed as a scattering length
    density (SLD).  This is the number density of the scatterers
    (per |Ang^3|) times their scattering lengths, scaled to
    |1e-6/Ang^2| (with |1/Ang^2| = $10^{5}$ fm/|Ang^3|).
    Following the convention of Sears (1992), we define sld as
    $\rho = \rho_{\rm re} - i \rho_{\rm im}$.

    .. math::

        \rho_{\rm re} (10^6 / Å^2) &= 10 N \mathrm{Re}(b_c) \\
        \rho_{\rm im} (10^6 / Å^2) &= -10 N \mathrm{Im}(b_c) \\
        \rho_{\rm inc} (10^6 / Å^2) &= 10 N b_i

    Similarly, the macroscopic scattering cross section for the sample includes
    number density:

    .. math::

        \Sigma_{\rm coh} (1/{\rm cm}) &= N \sigma_c \\
        \Sigma_{\rm inc} (1/{\rm cm}) &= N \sigma_i \\
        \Sigma_{\rm abs} (1/{\rm cm}) &= N \sigma_a \\
        \Sigma_{\rm s} (1/{\rm cm}) &= N \sigma_s

    The 1/e penetration depth $t_u$ represents the the depth into the sample
    at which the unscattered intensity is reduced by a factor of $e$:

    .. math::

        t_u (cm) = \left. 1 \right/ (\Sigma_{\rm s} + \Sigma_{\rm abs})

    The calculated penetration depth includes the effects of both
    absorption and incoherent scattering (which spreads the beam over the
    full $4\pi$ spherical surface, and so it looks like absorption with
    respect to the beam), as well as the coherent scattering from the sample.
    If you instead want to calculate the effective shielding of the sample,
    you should recalculate penetration depth with absorption only.

    Transmission rate can be computed from $e^{-d/t_u}$ for penetration
    depth $t_u$ and sample thickness $d$. This does not include many
    real world effects, such as single phonon scattering\ [#Mildner1998]_
    and forward scattering\ [#May1982]_, which result in measured
    transmission significantly different from the values predicted from
    nuclear properties alone.

    There is also a wavelength dependence for single phonon interactions which
    gives rise to significant inelastic scattering for lighter isotopes (H, D)
    and/or longer wavelengths (above 5 |Ang|). This factor is both
    temperature and material dependent and will not be included
    in the scattering calculations. In particular, penetration length and
    transmitted flux are going to be significantly overestimated.

    Including unit conversion with $\mu=10^{-6}$ the full scattering equations
    are:

    .. math::

        \rho_{\rm re}\,(\mu/Å^2) &= (N/Å^3)
            \, (\mathrm{Re}(b_c)\,{\rm fm})
            \, (10^{-5} Å/{\rm\,fm})
            \, (10^6\,\mu) \\
        \rho_{\rm im}\,(\mu/Å^2) &= (N/Å^3)
            \, (\sigma_a\,{\rm barn})
            \, (10^{-8}\,Å^2/{\rm barn}) / (2 \lambda\, Å)
            \, (10^6\,\mu) \\
            &= (N/Å^3)
            \, (-\mathrm{Im}(b_c)\,{\rm fm})
            \, (10^{-5} Å/{\rm\,fm})
            \, (10^6\,\mu) \\
        \rho_{\rm inc}\,(\mu/Å^2) &= (N/Å^3)
            \, \sqrt{(\sigma_i\, {\rm barn})/(4 \pi)
                \, (100\, {\rm fm}^2/{\rm barn})}
            \, (10^{-5}\, Å/{\rm fm})
            \, (10^6\, \mu) \\
        \Sigma_{\rm coh}\,(1/{\rm cm}) &= (N/Å^3)
            \, (\sigma_c\, {\rm barn})
            \, (10^{-8}\, Å^2/{\rm barn})
            \, (10^8\, Å/{\rm cm}) \\
        \Sigma_{\rm inc}\,(1/{\rm cm}) &= (N/Å^3)
            \,(\sigma_i\, {\rm barn})
            \, (10^{-8}\, Å^2/{\rm barn})
            \, (10^8\, Å/{\rm cm}) \\
        \Sigma_{\rm abs}\,(1/{\rm cm}) &= (N/Å^3)
            \,(\sigma_a\,{\rm barn})
            \, (10^{-8}\, Å^2/{\rm barn})
            \, (10^8\, Å/{\rm cm}) \\
        \Sigma_{\rm s}\,(1/{\rm cm}) &= (N/Å^3)
            \,(\sigma_s\,{\rm barn})
            \, (10^{-8}\, Å^2/{\rm barn})
            \, (10^8\, Å/{\rm cm}) \\
        t_u\,({\rm cm}) &= 1/(\Sigma_{\rm s}\, 1/{\rm cm}
            \,+\, \Sigma_{\rm abs}\, 1/{\rm cm})
    """
    from . import formulas

    compound = formulas.formula(
        compound, density=density, natural_density=natural_density, table=table)
    assert compound.density is not None, "scattering calculation needs density"
    #print("sld", compound, compound.density)
    if energy is not None:
        wavelength = neutron_wavelength(energy)
    # PAK: 1.5.3 wavelength now defaults to ABSORPTION_WAVELENGTH
    elif wavelength is None:
        wavelength = ABSORPTION_WAVELENGTH

    # Sum over the quantities (note: formulas can have fractional atoms)
    molar_mass, num_atoms = 0., 0.
    b_c, sigma_s = 0., 0.
    is_energy_dependent = False
    for element, quantity in compound.atoms.items():
        # TODO: use NaN rather than None
        if not element.neutron.has_sld():
            return None, None, None
        molar_mass += element.mass*quantity
        num_atoms += quantity
        # PAK 2021-04-05: allow energy dependent b_c, b''
        b_ck, sigma_sk = element.neutron.scattering_by_wavelength(wavelength)
        #print(f"{element=}; {b_ck=}; {sigma_sk=}")
        b_c += quantity * b_ck
        sigma_s += quantity * sigma_sk
        is_energy_dependent |= element.neutron.is_energy_dependent

    # If nothing to sum, return values for a vacuum.  This might be because
    # the material has no atoms or it might be because the density is zero.
    if molar_mass*compound.density == 0:
        return (0, 0, 0), (0, 0, 0), inf

    # Turn weighted sums into scattering factors
    b_c /= num_atoms
    sigma_s /= num_atoms

    # Compute number density (N/A^3)
    # volume A^3/N = ((1 g/mole) / (1 g/cm^3)) / (N/mole) * (10^8 A/cm)^3
    cell_volume = (molar_mass/compound.density)/avogadro_number*1e24
    number_density = num_atoms / cell_volume # N/A^3 = N/A^3

    return _calculate_scattering(number_density, wavelength, b_c, sigma_s)


def _calculate_scattering(number_density, wavelength, b_c, sigma_s):
    r"""
    :Parameters:
        *number_density* \: float | N/|Ang^3|
            Scatterers per unit volume.
        *wavelength* \: float(s) | |Ang|
            Neutron wavelength(s).
        *b_c* \: complex(s) | fm
            Complex bound coherent scattering length $b_c$.
        *sigma_c* \: float(s) | barn
            Total cross section.

    See neutron_scattering docstring for calculation details.

    Note: returns -sld_im for historical reasons.
    """
    #print("in scat", number_density, wavelength, b_c, sigma_s)
    # Compute SLD (1e-6/A^2). Extending Sears (1992) convention, b = b' - i b'',
    # returning sld = sld_re - i sld_im.
    sld = 10*number_density * b_c # 1e-6/A^2 = 1/A^3 1 fm 1e-5 A/fm 1e6 1e-6
    sld_re, sld_im = sld.real, abs(sld.imag)

    # PAK 2017-04-21: compute incoherent xs from total xs
    # PAK 2021-04-20: include imaginary b_c in coherent cross section
    # Compute coherent and incoherent cross sections (barn)
    sigma_c = _4PI_100 * abs(b_c)**2 # 1 barn = 1e-2 fm^2
    sigma_i = np.maximum(sigma_s - sigma_c, 0.)  # 1 barn = 1 barn

    # Compute incoherent scattering length from incoherent cross section (fm)
    b_i = sqrt(sigma_i / _4PI_100) # 1 fm = sqrt(1 barn * 1e-2 fm^2/barn)

    # Compute incoherent scattering length density (1e-6/A^2)
    sld_inc = number_density * b_i * 10 # 1e-6/A^2 = 1/A^3 1 fm 1e-5 A/fm 1e6 1e-6

    # Compute absorption cross section (barn)
    # Note: Sears (1992) uses b = b' - i b'', so use |Im(b_c)| for sigma_a.
    sigma_a = 2000 * abs(b_c.imag) * wavelength # 1 barn = 1 fm 1 A 1e5 A/fm 1e-2 barn/fm

    # print(f"σ_a {sigma_c:.3f} σ_i {sigma_i:.3f} σ_s {sigma_s:.3f} σ_a {sigma_a:.3f}")

    # Compute macroscopic scattering cross section per unit volume (1/cm)
    total_xs = number_density * sigma_s # 1/cm = 1/A^3 1 barn 1e-8 A^2/barn 1e8 A/cm
    coh_xs = number_density * sigma_c
    abs_xs = number_density * sigma_a
    inc_xs = number_density * sigma_i

    # Compute 1/e length (cm)
    penetration = 1/(abs_xs + total_xs) # 1 cm = 1 / (1/cm)

    return (sld_re, sld_im, sld_inc), (coh_xs, abs_xs, inc_xs), penetration


def neutron_sld(*args, **kw):
    r"""
    Computes neutron scattering length densities for molecules.

    :Parameters:
        *compound* \: Formula initializer
            Chemical formula
        *density* \: float | |g/cm^3|
            Mass density
        *natural_density* \: float | |g/cm^3|
            Mass density of formula with naturally occuring abundances
        *wavelength* \: float | |Ang|
            Neutron wavelength (default=1.798 |Ang|).
        *energy* \: float | meV
            Neutron energy.  If energy is specified then wavelength is ignored.
        *table* \: PeriodicTable
            Alternate table to use when parsing *compound*.
    :Returns:
        *sld* \: (float, float, float) | |1e-6/Ang^2|
            (*real*, -*imaginary*, *incoherent*) scattering length density.

    :Raises:
        *AssertionError* \: density is missing.

    Returns the scattering length density of the compound.
    See :func:`neutron_scattering` for details.
    """
    return neutron_scattering(*args, **kw)[0]

def neutron_sld_from_atoms(*args, **kw):
    r"""
    .. deprecated:: 0.91

        :func:`neutron_sld` accepts dictionaries of \{atom\: count\}.

    """
    return neutron_scattering(*args, **kw)[0]


def D2O_match(compound, **kw):
    """
    Find the D2O contrast match point for the compound.

    *wavelength* or *energy* select neutron wavelength or energy.

    Additional keyword arguments (*density*, *natural_density*, *name*, *table*)
    are passed to :func:`formulas.formula` when parsing the compound.

    Returns *D2O_fraction* and *SLD* at match point.

    See :func:`D2O_sld` for details on the calculation.

    Note that the resulting fraction is only meaningful in [0, 1]. Beyond
    100% you will need an additional constrast agent in the 100% D2O
    solvent to increase the SLD enough to match.
    """
    H2O_sld, D2O_sld, Hsld, Dsld = _D2O_slds(compound, **kw)
    # SLD(%Dsample + (1-%)Hsample) = SLD(%D2O + (1-%)H2O)
    # => %SLD(Dsample) + (1-%)SLD(Hsample) = %SLD(D2O) + (1-%)SLD(H2O)
    # => %(SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O))
    #      = SLD(H2O) - SLD(Hsample)
    # => % = 100*(SLD(H2O) - SLD(Hsample))
    #      / (SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O))
    D2O_fraction = (
        (H2O_sld[0] - Hsld[0]) / (Dsld[0] - Hsld[0] + H2O_sld[0] - D2O_sld[0]))

    match_point_sld = mix_values(Dsld, Hsld, D2O_fraction)
    return D2O_fraction, match_point_sld[0]

def D2O_sld(compound, volume_fraction=1., D2O_fraction=0., **kw):
    """
    Compute the neutron SLD for a D2O contrast solution.

    *compound* is a string or parsed formula object. Labile hydrogen should
    be marked as H[1] in the formula. These will be substituted according to
    %D2O in the solvent.

    The D2O contrast mixture is assumed to be made using pure H2O (with
    its natural H:D ratios) and pure D2O with no H present, so H[1] will be
    substituted alternately with H and D when computing mixture SLD.
    Solvent SLD is calculated using the density at 20 C.

    Only the coherent scattering crosssection will be matched. Incoherent
    and absorption crosssections are likely to be different for the compound
    and the solvent, especially due to the large incoherent crosssection for
    hydrogen.

    Note that incoherent scattering does not mix linearly, so the incoherent
    sld for the mixture will differ slightly from incoherent scattering
    computed returned from a compound with the same isotope ratios.

    *volume_fraction* is the portion by volume of solute in the solution.

    *D2O_fraction* is the portion by volume of D2O in the solvent.

    *wavelength* or *energy* to select neutron wavelength or energy.

    Additional keyword arguments (*density*, *natural_density*, *name*, *table*)
    are passed to :func:`formulas.formula` when parsing the compound.

    Returns (real, imag, incoh) SLD.
    """
    # TODO: fix incoherent scattering so it is consistent with compound
    # Need to compute sld from mixture rather than mixing parts
    H2O_sld, D2O_sld, Hsld, Dsld = _D2O_slds(compound, **kw)
    solvent_sld = mix_values(D2O_sld, H2O_sld, D2O_fraction)
    solute_sld = mix_values(Dsld, Hsld, D2O_fraction)
    solution_sld = mix_values(solute_sld, solvent_sld, volume_fraction)
    #print(D2O_fraction, volume_fraction)
    #print(compound, "solvent", solvent_sld)
    #print(compound, "solute", solute_sld)
    #print(compound, "solution", solution_sld)
    return solution_sld


def _D2O_slds(compound, **kw):
    from . import formulas

    # Water density at 20C; neutron wavelength doesn't matter.
    sld_args = dict(
        wavelength=kw.pop("wavelength", None),
        energy=kw.pop("energy", None),
        # Note: using get() rather than pop() for table since table can be a
        # parameter for formula and for neutron_sld (which calls formula)
        table=kw.get('table', None),
    )
    # TODO: use same table for solvent as solute?
    H2O_sld = neutron_sld("H2O@0.9982n", **sld_args)
    D2O_sld = neutron_sld("D2O@0.9982n", **sld_args)
    mol = formulas.formula(compound, **kw)
    # Be sure to pull H and H[1] from the table for the compound, otherwise
    # the elements may not match in the substitution.
    # TODO: include table in compound so parsed
    table = default_table(kw.get('table', None))
    labile_H, H, D = table.H[1], table.H, table.D
    Hsld = neutron_sld(mol.replace(labile_H, H), **sld_args)
    Dsld = neutron_sld(mol.replace(labile_H, D), **sld_args)

    return H2O_sld, D2O_sld, Hsld, Dsld


def mix_values(a, b, fraction):
    """
    Mix two tuples with floating point values according to fraction of a.
    """
    return tuple(aj*fraction + bj*(1-fraction) for aj, bj in zip(a, b))


def _sum_piece(wavelength, compound):
    """
    Helper for neutron_composite_sld which precomputes quantities of interest
    for material fragments in a composite formula.
    """
    # Sum over the quantities.
    molar_mass = num_atoms = 0
    b_c = sigma_s = 0
    for element, quantity in compound.atoms.items():
        molar_mass += element.mass*quantity
        num_atoms += quantity
        b_ck, sigma_sk = element.neutron.scattering_by_wavelength(wavelength)
        b_c += quantity * b_ck
        sigma_s += quantity * sigma_sk

    return num_atoms, molar_mass, b_c, sigma_s

# TODO: compute density from material densities if requested
# You ought to be able to use mixby="mass" or "volume" when creating the
# calculator, and ignore any density provided.
def neutron_composite_sld(materials, wavelength=ABSORPTION_WAVELENGTH):
    r"""
    Create a composite SLD calculator.

    :Parameters:
        *materials* \: [Formula]
            List of materials
        *wavelength* = 1.798: float OR [float] | |Ang|
            Probe wavelength(s).

    :Returns:
        *calculator* \: f(w, density=1) -> (*real*, -*imaginary*, *incoherent*)

    The composite calculator takes a vector of weights and returns the
    scattering length density of the composite.  This is useful for operations
    on large molecules, such as calculating a set of contrasts or fitting
    a material composition. Note that density must be provided for each set
    of material weights. The density on the individual materials is ignored.

    The returned slds will be vectors if the input wavelength is a vector
    and if any of the elements are energy dependent.

    Table lookups and partial sums and constants are precomputed so that
    the calculation consists of a few simple array operations regardless
    of the size of the material fragments.
    """
    # Input may be a scalar or a sequence. If it is a sequence, turn it into
    # an array before proceeding. If it is a scalar leave it as a scalar so
    # that float input returns float output.
    is_multi = not np.isscalar(wavelength)
    if is_multi:
        wavelength = np.asarray(wavelength)
    # Query all parts of the composition
    parts = [_sum_piece(wavelength, m) for m in materials]
    num_atoms_parts, molar_mass_parts, bc_parts, sigma_parts = [
        np.array(v) for v in zip(*parts)
    ]

    #for name, v in zip("N mass bc b'' total".split(), V): print(name, v)
    def _compute(weights, density=1):
        multiweights = weights[:, None] if is_multi else weights
        # Sum over the quantities
        molar_mass = np.sum(weights*molar_mass_parts)
        num_atoms = np.sum(weights*num_atoms_parts)
        b_c = np.sum(multiweights*bc_parts, axis=0)
        sigma_s = np.sum(multiweights*sigma_parts, axis=0)

        # If nothing to sum, return values for a vacuum.  This might be because
        # the material has no atoms or it might be because the density is zero.
        if molar_mass*density == 0:
            return 0, 0, 0

        # Compute number density (1/A^3)
        cell_volume = (molar_mass/density)/avogadro_number*1e24
        number_density = num_atoms / cell_volume
        #print("in compute", b_c, number_density)

        # Turn sums into scattering factors
        b_c /= num_atoms
        sigma_s /= num_atoms

        # TODO: duplicated from _calculate_scattering
        # Compute SLD (1e-6/A^2). Extending Sears (1992) convention, b = b' - i b'',
        # returning sld = sld_re - i sld_im.
        sld = 10*number_density * b_c # 1e-6/A^2 = 1/A^3 1 fm 1e-5 A/fm 1e6 1e-6
        sld_re, sld_im = sld.real, abs(sld.imag)

        # PAK 2017-04-21: compute incoherent xs from total xs
        # PAK 2021-04-20: include imaginary b_c in coherent cross section
        sigma_c = _4PI_100 * abs(b_c)**2 # 1 barn = 1e-2 fm^2
        sigma_i = np.maximum(sigma_s - sigma_c, 0.) # 1 barn = 1 barn
        b_i = sqrt(sigma_i / _4PI_100) # 1 fm = sqrt(1 barn * 1e-2 fm^2/barn)
        sld_inc = number_density * b_i * 10 # 1e-6/A^2 = 1/A^3 1 fm 1e-5 A/fm 1e6 1e-6

        return sld_re, sld_im, sld_inc

    return _compute


def sld_plot(table=None):
    r"""
    Plots SLD as a function of element number.

    :Parameters:
        *table* \: PeriodicTable
            The default periodictable unless a specific table has been requested.

    :Returns: None
    """
    from .plot import table_plot

    table = default_table(table)

    SLDs = dict((el, el.neutron.sld()[0])
                for el in table
                if el.neutron.has_sld())
    SLDs[table.D] = table.D.neutron.sld()[0]

    table_plot(SLDs, label='Scattering length density ($10^{-6}$ Nb)',
               title='Neutron SLD for elements in natural abundance')


# We are including the complete original table here in case somebody in
# future wants to extract uncertainties or other information.
#
# Z-Symbol-A
#   This is the atomic number, the symbol and the isotope.
#   If Z-Symbol only, the line represents an element with scattering determined
#   by the natural abundance of the isotopes in laboratory samples.  If there
#   is only one isotope, then there is no corresponding element definition.
# concentration/half-life
#   This is the natural abundance of the isotope expressed as a percentage, or
#   it is the half-life in years (number Y) or seconds (number S).
# spin I
#   For isotopes, the nuclear spin.
# b_c, bp, bm
#   Bound coherent scattering length in fm
#   b+/b- if present are spin dependent scattering for I+1/2 and I-1/2
#   respectively
# c
#   'E' if there is a strong energy dependency.
#   '+/-' if separate b+/b- values are available [PAK: doesn't seem true]
# coherent, incoherent, total
#   The coherent and incoherent scattering cross-sections in barns.
# absorption
#   The thermal absorption cross section in barns at 1.798 A; 25.30 meV; 2200 m/s
#
# Numbers in parenthesis represents uncertainty.
# Numbers followed by '*' are estimated.
# Numbers may be given as limit, e.g., <1.0e-6
#
# [Paul Kienzle]
# * Fix typos such as 70Zn b_c 6.9(1.0) => 6.0(1.0).
# * Update bound coherent scattering length for H-1, H-2, He-4, C-12,
#   O-16, O-17, O-18, Sn-119, Sm-154, Eu-153, Pb-207, Bi-209
# * Update total cross section for He, Kr, Xe
# * Usd 63-Eu-151 b_c from 84Mug1. This change is moot since this isotope
#   has energy dependent isotope coeffs.
# * Use calculated values for 4He coh and total, and natHe b_c, coh and total.

nsftable = """\
0-n-1,618 S,1/2,-37.0(6),0,-37.0(6),,43.01(2),,43.01(2),0
1-H,,,-3.7409(11),,,,1.7568(10),80.26(6),82.02(6),0.3326(7)
1-H-1,99.985,1/2,-3.7395(11),10.817(5),-47.420(14),+/-,1.7583(10),80.27(6),82.03(6),0.3326(7)
1-H-2,0.0149,1,6.6681(27),9.53(3),0.975(60),,5.592(7),2.05(3),7.64(3),0.000519(7)
1-H-3,12.26 Y,1/2,4.792(27),4.18(15),6.56(37),,2.89(3),0.14(4),3.03(5),<6.0E-6
2-He,,,3.0985(21),,,,1.2065(16),0,1.2065(16),0.00747(1)
2-He-3,0.013,1/2,5.74(7),4.7(5),8.8(1.4),E,4.42(10),1.6(4),6.0(4),5333.0(7.0)
2-He-4,99.987,0,3.0982(21),,,,1.2062(16),0,1.2062(16),0
3-Li,,,-1.90(3),,,,0.454(10),0.92(3),1.37(3),70.5(3)
3-Li-6,7.5,1,2.0(1),0.67(14),4.67(17),+/-,0.51(5),0.46(5),0.97(7),940.0(4.0)
3-Li-7,92.5,3/2,-2.22(2),-4.15(6),1.00(8),+/-,0.619(11),0.78(3),1.40(3),0.0454(3)
4-Be-9,100,3/2,7.79(1),,,,7.63(2),0.0018(9),7.63(2),0.0076(8)
5-B,,,5.30(4),,,,3.54(5),1.70(12),5.24(11),767.0(8.0)
5-B-10,19.4,3,-0.2(4),-4.2(4),5.2(4),,0.144(6),3.0(4),3.1(4),3835.0(9.0)
5-B-11,80.2,3/2,6.65(4),5.6(3),8.3(3),,5.56(7),0.21(7),5.77(10),0.0055(33)
6-C,,,6.6472(9),,,,5.551(2),0.001(4),5.551(3),0.00350(7)
6-C-12,98.89,0,6.6535(14),,,,5.559(3),0,5.559(3),0.00353(7)
6-C-13,1.11,1/2,6.19(9),5.6(5),6.2(5),+/-,4.81(14),0.034(11),4.84(14),0.00137(4)
7-N,,,9.36(2),,,,11.01(5),0.50(12),11.51(11),1.90(3)
7-N-14,99.635,1,9.37(2),10.7(2),6.2(3),,11.03(5),0.50(12),11.53(11),1.91(3)
7-N-15,0.365,1/2,6.44(3),6.77(10),6.21(10),,5.21(5),0.00005(10),5.21(5),0.000024(8)
8-O,,,5.8037(29),,,,4.232(6),0.000(8),4.232(6),0.00019(2)
8-O-16,99.75,0,5.805(5),,,,4.232(6),0,4.232(6),0.00010(2)
8-O-17,0.039,5/2,5.867(4),5.52(20),5.17(20),,4.20(22),0.004(3),4.20(22),0.236(10)
8-O-18,0.208,0,6.009(5),,,,4.29(10),0,4.29(10),0.00016(1)
9-F-19,100,1/2,5.654(12),5.632(10),5.767(10),+/-,4.017(14),0.0008(2),4.018(14),0.0096(5)
10-Ne,,,4.566(6),,,,2.620(7),0.008(9),2.628(6),0.039(4)
10-Ne-20,90.5,0,4.631(6),,,,2.695(7),0,2.695(7),0.036(4)
10-Ne-21,0.27,3/2,6.66(19),,,,5.6(3),0.05(2),5.7(3),0.67(11)
10-Ne-22,9.2,0,3.87(1),,,,1.88(1),0,1.88(1),0.046(6)
11-Na-23,100,3/2,3.63(2),6.42(4),-1.00(6),+/-,1.66(2),1.62(3),3.28(4),0.530(5)
12-Mg,,,5.375(4),,,,3.631(5),0.08(6),3.71(4),0.063(3)
12-Mg-24,78.99,0,5.49(18),,,,4.03(4),0,4.03(4),0.050(5)
12-Mg-25,10,5/2,3.62(14),4.73(30),1.76(20),+/-,1.65(13),0.28(4),1.93(14),0.19(3)
12-Mg-26,11,0,4.89(15),,,,3.00(18),0,3.00(18),0.0382(8)
13-Al-27,100,5/2,3.449(5),3.67(2),3.15(2),,1.495(4),0.0082(6),1.503(4),0.231(3)
14-Si,,,4.15071(22),,,,2.1633(10),0.004(8),2.167(8),0.171(3)
14-Si-28,92.2,0,4.106(6),,,,2.120(6),0,2.120(6),0.177(3)
14-Si-29,4.7,1/2,4.7(1),4.50(15),4.7(4),+/-,2.78(12),0.001(2),2.78(12),0.101(14)
14-Si-30,3.1,0,4.58(8),,,,2.64(9),0,2.64(9),0.107(2)
15-P-31,100,1/2,5.13(1),,,+/-,3.307(13),0.005(10),3.312(16),0.172(6)
16-S,,,2.847(1),,,,1.0186(7),0.007(5),1.026(5),0.53(1)
16-S-32,95,0,2.804(2),,,,0.9880(14),0,0.9880(14),0.54(4)
16-S-33,0.74,3/2,4.74(19),,,+/-,2.8(2),0.3(6),3.1(6),0.54(4)
16-S-34,4.2,0,3.48(3),,,,1.52(3),0,1.52(3),0.227(5)
16-S-36,0.02,0,3.0(1.0)*,,,,1.1(8),0,1.1(8),0.15(3)
17-Cl,,,9.5792(8),,,,11.528(2),5.3(5),16.8(5),33.5(3)
17-Cl-35,75.77,3/2,11.70(9),16.3(2),4.0(3),+/-,17.06(6),4.7(6),21.8(6),44.1(4)
17-Cl-37,24.23,3/2,3.08(6),3.10(7),3.05(7),+/-,1.19(5),0.001(3),1.19(5),0.433(6)
18-Ar,,,1.909(6),,,,0.458(3),0.225(5),0.683(4),0.675(9)
18-Ar-36,0.34,0,24.9(7),,,,77.9(4),0,77.9(4),5.2(5)
18-Ar-38,0.07,0,3.5(3.5),,,,1.5(3.1),0,1.5(3.1),0.8(5)
18-Ar-40,99.59,0,1.7,,,,0.421(3),0,0.421(3),0.660(9)
19-K,,,3.67(2),,,,1.69(2),0.27(11),1.96(11),2.1(1)
19-K-39,93.3,3/2,3.79(2),5.15,1.51,+/-,1.76(2),0.25(11),2.01(11),2.1(1)
19-K-40,0.012,4,3.1(1.0)*,,,,1.1(6),0.5(5)*,1.6(9),35.0(8.0)
19-K-41,6.7,3/2,2.69(8),,,,0.91(5),0.3(6),1.2(6),1.46(3)
20-Ca,,,4.70(2),,,,2.78(2),0.05(3),2.83(2),0.43(2)
20-Ca-40,96.94,0,4.78(5),,,,2.90(2),0,2.90(2),0.41(2)
20-Ca-42,0.64,0,3.36(10),,,,1.42(8),0,1.42(8),0.68(7)
20-Ca-43,0.13,7/2,-1.56(9),,,,0.31(4),0.5(5),0.8(5),6.2(6)
20-Ca-44,2.13,0,1.42(6),,,,0.25(2),0,0.25(2),0.88(5)
20-Ca-46,0.003,0,3.55(21),,,,1.6(2),0,1.6(2),0.74(7)
20-Ca-48,0.18,0,0.39(9),,,,0.019(9),0,0.019(9),1.09(14)
21-Sc-45,100,7/2,12.1(1),6.91(22),18.99(28),+/-,19.0(3),4.5(3),23.5(6),27.5(2)
22-Ti,,,-3.370(13),,,,1.485(2),2.87(3),4.35(3),6.09(13)
22-Ti-46,8,0,4.72(5),,,,3.05(7),0,3.05(7),0.59(18)
22-Ti-47,7.5,5/2,3.53(7),0.46(23),7.64(13),,1.66(11),1.5(2),3.2(2),1.7(2)
22-Ti-48,73.7,0,-5.86(2),,,,4.65(3),0,4.65(3),7.84(25)
22-Ti-49,5.5,7/2,0.98(5),2.6(3),-1.2(4),,0.14(1),3.3(3),3.4(3),2.2(3)
22-Ti-50,5.3,0,5.88(10),,,,4.80(12),0,4.80(12),0.179(3)
23-V,,,-0.443(14),,,,0.01838(12),5.08(6),5.10(6),5.08(4)
23-V-50,0.25,6,7.6(6)*,,,,7.3(1.1),0.5(5)*,7.8(1.0),60.0(40.0)
23-V-51,99.75,7/2,-0.402(2),4.93(25),-7.58(28),+/-,0.0203(2),5.07(6),5.09(6),4.9(1)
24-Cr,,,3.635(7),,,,1.660(6),1.83(2),3.49(2),3.05(6)
24-Cr-50,4.35,0,-4.50(5),,,,2.54(6),0,2.54(6),15.8(2)
24-Cr-52,83.8,0,4.914(15),,,,3.042(12),0,3.042(12),0.76(6)
24-Cr-53,9.59,3/2,-4.20(3),1.16(10),-13.0(2),,2.22(3),5.93(17),8.15(17),18.1(1.5)
24-Cr-54,2.36,0,4.55(10),,,,2.60(11),0,2.60(11),0.36(4)
25-Mn-55,100,5/2,-3.750(18),-4.93(46),-1.46(33),,1.75(2),0.40(2),2.15(3),13.3(2)
26-Fe,,,9.45(2),,,,11.22(5),0.40(11),11.62(10),2.56(3)
26-Fe-54,5.8,0,4.2(1),,,,2.2(1),0,2.2(1),2.25(18)
26-Fe-56,91.7,0,10.1(2),,,,12.42(7),0,12.42(7),2.59(14)
26-Fe-57,2.19,1/2,2.3(1),,,,0.66(6),0.3(3)*,1.0(3),2.48(30)
26-Fe-58,0.28,0,15(7),,,,28.0(26.0),0,28.0(26.0),1.28(5)
27-Co-59,100,7/2,2.49(2),-9.21(10),3.58(10),+/-,0.779(13),4.8(3),5.6(3),37.18(6)
28-Ni,,,10.3(1),,,,13.3(3),5.2(4),18.5(3),4.49(16)
28-Ni-58,67.88,0,14.4(1),,,,26.1(4),0,26.1(4),4.6(3)
28-Ni-60,26.23,0,2.8(1),,,,0.99(7),0,0.99(7),2.9(2)
28-Ni-61,1.19,3/2,7.60(6),,,,7.26(11),1.9(3),9.2(3),2.5(8)
28-Ni-62,3.66,0,-8.7(2),,,,9.5(4),0,9.5(4),14.5(3)
28-Ni-64,1.08,0,-0.37(7),,,,0.017(7),0,0.017(7),1.52(3)
29-Cu,,,7.718(4),,,,7.485(8),0.55(3),8.03(3),3.78(2)
29-Cu-63,69.1,3/2,6.477(13),,,+/-,5.2(2),0.006(1),5.2(2),4.50(2)
29-Cu-65,30.9,3/2,10.204(20),,,+/-,14.1(5),0.40(4),14.5(5),2.17(3)
30-Zn,,,5.680(5),,,,4.054(7),0.077(7),4.131(10),1.11(2)
30-Zn-64,48.9,0,5.23(4),,,,3.42(5),0,3.42(5),0.93(9)
30-Zn-66,27.8,0,5.98(5),,,,4.48(8),0,4.48(8),0.62(6)
30-Zn-67,4.1,5/2,7.58(8),5.8(5),10.1(7),+/-,7.18(15),0.28(3),7.46(15),6.8(8)
30-Zn-68,18.6,0,6.04(3),,,,4.57(5),0,4.57(5),1.1(1)
30-Zn-70,0.62,0,6.0(1.0)*,,,,4.5(1.5),0,4.5(1.5),0.092(5)
31-Ga,,,7.288(2),,,,6.675(4),0.16(3),6.83(3),2.75(3)
31-Ga-69,60,3/2,8.043(16),6.3(2),10.5(4),+/-,7.80(4),0.091(11),7.89(4),2.18(5)
31-Ga-71,40,3/2,6.170(11),5.5(6),7.8(1),+/-,5.15(5),0.084(8),5.23(5),3.61(10)
32-Ge,,,8.185(20),,,,8.42(4),0.18(7),8.60(6),2.20(4)
32-Ge-70,20.7,0,10.0(1),,,,12.6(3),0,12.6(3),3.0(2)
32-Ge-72,27.5,0,8.51(10),,,,9.1(2),0,9.1(2),0.8(2)
32-Ge-73,7.7,9/2,5.02(4),8.1(4),1.2(4),,3.17(5),1.5(3),4.7(3),15.1(4)
32-Ge-74,36.4,0,7.58(10),,,,7.2(2),0,7.2(2),0.4(2)
32-Ge-76,7.7,0,8.2(1.5),,,,8.0(3.0),0,8.0(3.0),0.16(2)
33-As-75,100,3/2,6.58(1),6.04(5),7.47(8),+/-,5.44(2),0.060(10),5.50(2),4.5(1)
34-Se,,,7.970(9),,,,7.98(2),0.32(6),8.30(6),11.7(2)
34-Se-74,0.9,0,0.8(3.0),,,,0.1(6),0,0.1(6),51.8(1.2)
34-Se-76,9,0,12.2(1),,,,18.7(3),0,18.7(3),85.0(7.0)
34-Se-77,7.5,0,8.25(8),,,,8.6(2),0.05(25),8.65(16),42.0(4.0)
34-Se-78,23.5,0,8.24(9),,,,8.5(2),0,8.5(2),0.43(2)
34-Se-80,50,0,7.48(3),,,,7.03(6),0,7.03(6),0.61(5)
34-Se-82,8.84,0,6.34(8),,,,5.05(13),0,5.05(13),0.044(3)
35-Br,,,6.79(2),,,,5.80(3),0.10(9),5.90(9),6.9(2)
35-Br-79,50.49,3/2,6.79(7),,,+/-,5.81(2),0.15(6),5.96(13),11.0(7)
35-Br-81,49.31,3/2,6.78(7),,,+/-,5.79(12),0.05(2),5.84(12),2.7(2)
36-Kr,,,7.81(2),,,,7.67(4),0.01(14),7.685(26),25.0(1.0)
36-Kr-78,0.35,0,,,,,,0,,6.4(9)
36-Kr-80,2.5,0,,,,,,0,,11.8(5)
36-Kr-82,11.6,0,,,,,,0,,29.0(20.0)
36-Kr-83,11.5,9/2,,,,,,,,185.0(30.0)
36-Kr-84,57,0,,,,,,0,6.6,0.113(15)
36-Kr-86,17.3,0,8.07(26),,,,8.2(4),0,8.2(4),0.003(2)
37-Rb,,,7.08(2),,,,6.32(4),0.5(4),6.8(4),0.38(1)
37-Rb-85,72.17,5/2,7.07(10),,,,6.2(2),0.5(5)*,6.7(5),0.48(1)
37-Rb-87,27.83,3/2,7.27(12),,,,6.6(2),0.5(5)*,7.1(5),0.12(3)
38-Sr,,,7.02(2),,,,6.19(4),0.06(11),6.25(10),1.28(6)
38-Sr-84,0.56,0,5.0(2.0),,,,6.0(2.0),0,6.0(2.0),0.87(7)
38-Sr-86,9.9,0,5.68(5),,,,4.04(7),0,4.04(7),1.04(7)
38-Sr-87,7,9/2,7.41(7),,,,6.88(13),0.5(5)*,7.4(5),16.0(3.0)
38-Sr-88,82.6,0,7.16(6),,,,6.42(11),0,6.42(11),0.058(4)
39-Y-89,100,1/2,7.75(2),8.4(2),5.8(5),+/-,7.55(4),0.15(8),7.70(9),1.28(2)
40-Zr,,,7.16(3),,,,6.44(5),0.02(15),6.46(14),0.185(3)
40-Zr-90,51.48,0,6.5(1),,,,5.1(2),0,5.1(2),0.011(5)
40-Zr-91,11.23,5/2,8.8(1),7.9(2),10.1(2),+/-,9.5(2),0.15(4),9.7(2),1.17(10)
40-Zr-92,17.11,0,7.5(2),,,,6.9(4),0,6.9(4),0.22(6)
40-Zr-94,17.4,0,8.3(2),,,,8.4(4),0,8.4(4),0.0499(24)
40-Zr-96,2.8,0,5.5(1),,,,3.8(1),0,3.8(1),0.0229(10)
41-Nb-93,100,9/2,7.054(3),7.06(4),7.35(4),+/-,6.253(5),0.0024(3),6.255(5),1.15(6)
42-Mo,,,6.715(20),,,,5.67(3),0.04(5),5.71(4),2.48(4)
42-Mo-92,15.48,0,6.93(8),,,,6.00(14),0,6.00(14),0.019(2)
42-Mo-94,9.1,0,6.82(7),,,,5.81(12),0,5.81(12),0.015(2)
42-Mo-95,15.72,5/2,6.93(7),,,,6.00(10),0.5(5)*,6.5(5),13.1(3)
42-Mo-96,16.53,0,6.22(6),,,,4.83(9),0,4.83(9),0.5(2)
42-Mo-97,9.5,5/2,7.26(8),,,,6.59(15),0.5(5)*,7.1(5),2.5(2)
42-Mo-98,23.78,0,6.60(7),,,,5.44(12),0,5.44(12),0.127(6)
42-Mo-100,9.6,0,6.75(7),,,,5.69(12),0,5.69(12),0.4(2)
43-Tc-99,210000 Y,9/2,6.8(3),,,,5.8(5),0.5(5)*,6.3(7),20.0(1.0)
44-Ru,,,7.02(2),,,,6.21(5),0.4(1),6.6(1),2.56(13)
44-Ru-96,5.8,0,,,,,,0,,0.28(2)
44-Ru-98,1.9,0,,,,,,0,,<8.0
44-Ru-99,12.7,5/2,,,,,,,,6.9(1.0)
44-Ru-100,12.6,0,,,,,,0,,4.8(6)
44-Ru-101,17.07,5/2,,,,,,,,3.3(9)
44-Ru-102,31.61,0,,,,,,0,,1.17(7)
44-Ru-104,18.58,0,,,,,,0,,0.31(2)
45-Rh-103,100,1/2,5.90(4),8.15(6),6.74(6),,4.34(6),0.3(3)*,4.6(3),144.8(7)
46-Pd,,,5.91(6),,,,4.39(9),0.093(9),4.48(9),6.9(4)
46-Pd-102,1,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),3.4(3)
46-Pd-104,11,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),0.6(3)
46-Pd-105,22.33,5/2,5.5(3),,,+/-,3.8(4),0.8(1.0),4.6(1.1),20.0(3.0)
46-Pd-106,27.33,0,6.4(4),,,,5.1(6),0,5.1(6),0.304(29)
46-Pd-108,26.71,0,4.1(3),,,,2.1(3),0,2.1(3),8.5(5)
46-Pd-110,11.8,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),0.226(31)
47-Ag,,,5.922(7),,,,4.407(10),0.58(3),4.99(3),63.3(4)
47-Ag-107,51.8,1/2,7.555(11),8.14(9),5.8(3),+/-,7.17(2),0.13(3),7.30(4),37.6(1.2)
47-Ag-109,48.2,1/2,4.165(11),3.24(8),6.9(2),+/-,2.18(1),0.32(5),2.50(5),91.0(1.0)
48-Cd,,,4.83(5),,,E,3.04(6),3.46(13),6.50(12),2520.0(50.0)
48-Cd-106,1.2,0,5.0(2.0)*,,,,3.1(2.5),0,3.1(2.5),1.0(2.0)
48-Cd-108,0.9,0,5.31(24),,,,3.7(1),0,3.7(1),1.1(3)
48-Cd-110,12.39,0,5.78(8),,,,4.4(1),0,4.4(1),11.0(1.0)
48-Cd-111,12.75,1/2,6.47(8),,,,5.3(2),0.3(3)*,5.6(4),24.0(5.0)
48-Cd-112,24.07,0,6.34(6),,,,5.1(2),0,5.1(2),2.2(5)
48-Cd-113,12.36,1/2,-8.0(1),,,E,12.1(4),0.3(3)*,12.4(5),20600.0(400.0)
48-Cd-114,28.86,0,7.48(5),,,,7.1(2),0,7.1(2),0.34(2)
48-Cd-116,7.58,0,6.26(9),,,,5.0(2),0,5.0(2),0.075(13)
49-In,,,4.065(20),,,,2.08(2),0.54(11),2.62(11),193.8(1.5)
49-In-113,4.28,9/2,5.39(6),,,,3.65(8),0.000037(5),3.65(8),12.0(1.1)
49-In-115,95.72,9/2,4.00(3),2.1(1),6.4(4),,2.02(2),0.55(11),2.57(11),202.0(2.0)
50-Sn,,,6.2239(13),,,,4.871(3),0.022(5),4.892(6),0.626(9)
50-Sn-112,1,0,6.0(1.0)*,,,,4.5(1.5),0,4.5(1.5),1.00(11)
50-Sn-114,0.66,0,6.0(3),,,,4.8(5),0,4.8(5),0.114(30)
50-Sn-115,0.35,1/2,6.0(1.0)*,,,,4.5(1.5),0.3(3)*,4.8(1.5),30.0(7.0)
50-Sn-116,14.3,0,6.10(1),,,,4.42(7),0,4.42(7),0.14(3)
50-Sn-117,7.61,1/2,6.59(8),0.22(10),-0.23(10),,5.28(8),0.3(3)*,5.6(3),2.3(5)
50-Sn-118,24.03,0,6.23(4),,,,4.63(8),0,4.63(8),0.22(5)
50-Sn-119,8.58,1/2,6.28(3),0.14(10),0.0(1),,4.71(8),0.3(3)*,5.0(3),2.2(5)
50-Sn-120,32.86,0,6.67(4),,,,5.29(8),0,5.29(8),0.14(3)
50-Sn-122,4.72,0,5.93(3),,,,4.14(7),0,4.14(7),0.18(2)
50-Sn-124,5.94,0,6.15(3),,,,4.48(8),0,4.48(8),0.133(5)
51-Sb,,,5.57(3),,,,3.90(4),0.00(7),3.90(6),4.91(5)
51-Sb-121,57.25,5/2,5.71(6),5.7(2),5.8(2),,4.10(9),0.0003(19),4.10(19),5.75(12)
51-Sb-123,42.75,7/2,5.38(7),5.2(2),5.4(2),,3.64(9),0.001(4),3.64(9),3.8(2)
52-Te,,,5.68(2),,,,4.23(4),0.09(6),4.32(5),4.7(1)
52-Te-120,0.09,0,5.3(5),,,,3.5(7),0,3.5(7),2.3(3)
52-Te-122,2.4,0,3.8(2),,,,1.8(2),0,1.8(2),3.4(5)
52-Te-123,0.87,1/2,-0.05(25),-1.2(2),3.5(2),,0.002(3),0.52(5),0.52(5),418.0(30.0)
52-Te-124,4.61,0,7.95(10),,,,8.0(2),0,8.0(2),6.8(1.3)
52-Te-125,6.99,1/2,5.01(8),4.9(2),5.5(2),,3.17(10),0.008(8),3.18(10),1.55(16)
52-Te-126,18.71,0,5.55(7),,,,3.88(10),0,3.88(10),1.04(15)
52-Te-128,31.79,0,5.88(8),,,,4.36(10),0,4.36(10),0.215(8)
52-Te-130,34.48,0,6.01(7),,,,4.55(11),0,4.55(11),0.29(6)
53-I-127,100,5/2,5.28(2),6.6(2),3.4(2),,3.50(3),0.31(6),3.81(7),6.15(6)
54-Xe,,,4.69(4),,,,3.04(4),0,4.344(17),23.9(1.2)
54-Xe-124,0.1,0,,,,,,0,,165.0(20.0)
54-Xe-126,0.09,0,,,,,,0,,3.5(8)
54-Xe-128,1.9,0,,,,,,0,,<8.0
54-Xe-129,26.14,1/2,,,,,,,,21.0(5.0)
54-Xe-130,3.3,0,,,,,,0,,<26.0
54-Xe-131,21.18,3/2,,,,,,,,85.0(10.0)
54-Xe-132,26.89,0,,,,,,0,,0.45(6)
54-Xe-134,10.4,0,,,,,,0,,0.265(20)
54-Xe-136,8.9,0,,,,,,0,,0.26(2)
55-Cs-133,100,7/2,5.42(2),,,+/-,3.69(15),0.21(5),3.90(6),29.0(1.5)
56-Ba,,,5.07(3),,,,3.23(4),0.15(11),3.38(10),1.1(1)
56-Ba-130,0.1,0,-3.6(6),,,,1.6(5),0,1.6(5),30.0(5.0)
56-Ba-132,0.09,0,7.8(3),,,,7.6(6),0,7.6(6),7.0(8)
56-Ba-134,2.4,0,5.7(1),,,,4.08(14),0,4.08(14),2.0(1.6)
56-Ba-135,6.59,3/2,4.66(10),,,,2.74(12),0.5(5)*,3.2(5),5.8(9)
56-Ba-136,7.81,0,4.90(8),,,,3.03(10),0,3.03(10),0.68(17)
56-Ba-137,11.32,3/2,6.82(10),,,,5.86(17),0.5(5)*,6.4(5),3.6(2)
56-Ba-138,71.66,0,4.83(8),,,,2.94(10),0,2.94(10),0.27(14)
57-La,,,8.24(4),,,,8.53(8),1.13(19),9.66(17),8.97(2)
57-La-138,0.09,5,8.0(2.0)*,,,,8.0(4.0),0.5(5)*,8.5(4.0),57.0(6.0)
57-La-139,99.91,7/2,8.24(4),11.4(3),4.5(4),+/-,8.53(8),1.13(15),9.66(17),8.93(4)
58-Ce,,,4.84(2),,,,2.94(2),0.00(10),2.94(10),0.63(4)
58-Ce-136,0.19,0,5.76(9),,,,4.23(13),0,4.23(13),7.3(1.5)
58-Ce-138,0.26,0,6.65(9),,,,5.64(15),0,5.64(15),1.1(3)
58-Ce-140,88.48,0,4.81(9),,,,2.94(11),0,2.94(11),0.57(4)
58-Ce-142,11.07,0,4.72(9),,,,2.84(11),0,2.84(11),0.95(5)
59-Pr-141,100,5/2,4.58(5),,,+/-,2.64(6),0.015(3),2.66(6),11.5(3)
60-Nd,,,7.69(5),,,,7.43(19),9.2(8),16.6(8),50.5(1.2)
60-Nd-142,27.11,0,7.7(3),,,,7.5(6),0,7.5(6),18.7(7)
60-Nd-143,12.17,7/2,14.0(2.0)*,,,,25.0(7.0),55.0(7.0),80.0(2.0),337.0(10.0)
60-Nd-144,23.85,0,2.8(3),,,,1.0(2),0,1.0(2),3.6(3)
60-Nd-145,8.5,7/2,14.0(2.0)*,,,,25.0(7.0),5.0(5.0)*,30.0(9.0),42.0(2.0)
60-Nd-146,17.22,0,8.7(2),,,,9.5(4),0,9.5(4),1.4(1)
60-Nd-148,5.7,0,5.7(3),,,,4.1(4),0,4.1(4),2.5(2)
60-Nd-150,5.6,0,5.28(20),,,,3.5(3),0,3.5(3),1.2(2)
61-Pm-147,2.62 Y,7/2,12.6(4),,,,20.0(1.3),1.3(2.0),21.3(1.5),168.4(3.5)
62-Sm,,,0.00(5),,,E,0.422(9),39.0(3.0),39.4(3.0),5922.0(56.0)
62-Sm-144,3.1,0,-3.0(4.0)*,,,,1.0(3.0),0,1.0(3.0),0.7(3)
62-Sm-147,15,7/2,14.0(3.0),,,,25.0(11.0),14.0(19.0),39.0(16.0),57.0(3.0)
62-Sm-148,11.2,0,-3.0(4.0)*,,,,1.0(3.0),0,1.0(3.0),2.4(6)
62-Sm-149,13.8,7/2,18.7(28),,,E,63.5(6),137.0(5.0),200.0(5.0),42080.0(400.0)
62-Sm-150,7.4,0,14.0(3.0),,,,25.0(11.0),0,25.0(11.0),104.0(4.0)
62-Sm-152,26.7,0,-5.0(6),,,,3.1(8),0,3.1(8),206.0(6.0)
62-Sm-154,22.8,0,8.97(6),,,,11.0(2.0),0,11.0(2.0),8.4(5)
63-Eu,,,5.3(3),,,E,6.57(4),2.5(4),9.2(4),4530.0(40.0)
63-Eu-151,47.8,5/2,6.92(15),,,E,5.5(2),3.1(4),8.6(4),9100.0(100.0)
63-Eu-153,52.8,5/2,8.85(3),,,,8.5(2),1.3(7),9.8(7),312.0(7.0)
64-Gd,,,9.5(2),,,E,29.3(8),151.0(2.0),180.0(2.0),49700.0(125.0)
64-Gd-152,0.2,0,10.0(3.0)*,,,,13.0(8.0),0,13.0(8.0),735.0(20.0)
64-Gd-154,2.2,0,10.0(3.0)*,,,,13.0(8.0),0,13.0(8.0),85.0(12.0)
64-Gd-155,14.9,3/2,13.8(3),,,E,40.8(4),25.0(6.0),66.0(6.0),61100.0(400.0)
64-Gd-156,20.6,0,6.3(4),,,,5.0(6),0,5.0(6),1.5(1.2)
64-Gd-157,15.7,3/2,4.0(2.0),,,E,650.0(4.0),394.0(7.0),1044.0(8.0),259000.0(700.0)
64-Gd-158,24.7,0,9.0(2.0),,,,10.0(5.0),0,10.0(5.0),2.2(2)
64-Gd-160,21.7,0,9.15(5),,,,10.52(11),0,10.52(11),0.77(2)
65-Tb-159,100,3/2,7.34(2),6.8(2),8.1(2),+/-,6.84(6),0.004(3),6.84(6),23.4(4)
66-Dy,,,16.9(3),,,,35.9(8),54.4(1.2),90.3(9),994.0(13.0)
66-Dy-156,0.06,0,6.1(5),,,,4.7(8),0,4.7(8),33.0(3.0)
66-Dy-158,0.1,0,6.0(4.0)*,,,,5.0(6.0),0,5.(6.),43.0(6.0)
66-Dy-160,2.3,0,6.7(4),,,,5.6(7),0,5.6(7),56.0(5.0)
66-Dy-161,18.9,5/2,10.3(4),,,,13.3(1.0),3.0(1.0),16.0(1.0),600.0(25.0)
66-Dy-162,25.5,0,-1.4(5),,,,0.25(18),0,0.25(18),194.0(10.0)
66-Dy-163,24.9,5/2,5.0(4),6.1(5),3.5(5),,3.1(5),0.21(19),3.3(5),124.0(7.0)
66-Dy-164,28.2,0,49.4(5),,,,307.0(3.0),0,307.0(3.0),2840.0(40.0)
67-Ho-165,100,7/2,8.44(3),6.9(2),10.3(2),+/-,8.06(8),0.36(3),8.42(16),64.7(1.2)
68-Er,,,7.79(2),,,,7.63(4),1.1(3),8.7(3),159.0(4.0)
68-Er-162,0.14,0,9.01(11),,,,9.7(4),0,9.7(4),19.0(2.0)
68-Er-164,1.6,0,7.95(14),,,,8.4(4),0,8.4(4),13.0(2.0)
68-Er-166,33.4,0,10.51(19),,,,14.1(5),0,14.1(5),19.6(1.5)
68-Er-167,22.9,7/2,3.06(5),5.3(3),0.0(3),,1.1(2),0.13(6),1.2(2),659.0(16.0)
68-Er-168,27,0,7.43(8),,,,6.9(7),0,6.9(7),2.74(8)
68-Er-170,15,0,9.61(6),,,,11.6(1.2),0,11.6(1.2),5.8(3)
69-Tm-169,100,1/2,7.07(3),,,+/-,6.28(5),0.10(7),6.38(9),100.0(2.0)
70-Yb,,,12.41(3),,,,19.42(9),4.0(2),23.4(2),34.8(8)
70-Yb-168,0.14,0,-4.07(2),,,E,2.13(2),0,2.13(2),2230.0(40.0)
70-Yb-170,3,0,6.8(1),,,,5.8(2),0,5.8(2),11.4(1.0)
70-Yb-171,14.3,1/2,9.7(1),6.5(2),19.4(4),,11.7(2),3.9(2),15.6(3),48.6(2.5)
70-Yb-172,21.9,0,9.5(1),,,,11.2(2),0,11.2(2),0.8(4)
70-Yb-173,16.3,5/2,9.56(10),2.5(2),13.3(3),,11.5(2),3.5,15,17.1(1.3)
70-Yb-174,31.8,0,19.2(1),,,,46.8(5),0,46.8(5),69.4(5.0)
70-Yb-176,12.7,0,8.7(1),,,,9.6(2),0,9.6(2),2.85(5)
71-Lu,,,7.21(3),,,,6.53(5),0.7(4),7.2(4),74.0(2.0)
71-Lu-175,97.4,7/2,7.28(9),,,,6.59(5),0.6(4),7.2(4),21.0(3.0)
71-Lu-176,2.6,7,6.1(2),,,,4.7(2),1.2(3),5.9,2065.(35.)
72-Hf,,,7.77(14),,,,7.6(3),2.6(5),10.2(4),104.1(5)
72-Hf-174,0.184,0,10.9(1.1),,,,15.0(3.0),0,15.0(3.0),561.0(35.0)
72-Hf-176,5.2,0,6.61(18),,,,5.5(3),0,5.5(3),23.5(3.1)
72-Hf-177,18.5,0,0.8(1.0)*,,,,0.1(2),0.1(3),0.2(2),373.0(10.0)
72-Hf-178,27.2,0,5.9(2),,,,4.4(3),0,4.4(3),84.0(4.0)
72-Hf-179,13.8,9/2,7.46(16),,,,7.0(3),0.14(2),7.1(3),41.0(3.0)
72-Hf-180,35.1,0,13.2(3),,,,21.9(1.0),0,21.9(1.0),13.04(7)
73-Ta,,,6.91(7),,,,6.00(12),0.01(17),6.01(12),20.6(5)
73-Ta-180,0.012,9,7.0(2.0)*,,,,6.2(3.5),0.5(5)*,7.0(4.0),563.0(60.0)
73-Ta-181,99.98,7/2,6.91(7),,,+/-,6.00(12),0.011(2),6.01(12),20.5(5)
74-W,,,4.755(18),,,,2.97(2),1.63(6),4.60(6),18.3(2)
74-W-180,0.13,0,5.0(3.0)*,,,,3.0(4.0),0,3.0(4.0),30.0(20.0)
74-W-182,26.3,1/2,7.04(4),,,,6.10(7),0,6.10(7),20.7(5)
74-W-183,14.3,1/2,6.59(4),6.3(4),7.0(4),,5.36(7),0.3(3)*,5.7(3),10.1(3)
74-W-184,30.7,0,7.55(6),,,,7.03(11),0,7.03(11),1.7(1)
74-W-186,28.6,0,-0.73(4),,,,0.065(7),0,0.065(7),37.9(6)
75-Re,,,9.2(2),,,,10.6(5),0.9(6),11.5(3),89.7(1.0)
75-Re-185,37.5,5/2,9.0(3),,,,10.2(7),0.5(9),10.7(6),112.0(2.0)
75-Re-187,62.5,5/2,9.3(3),,,,10.9(7),1.0(6),11.9(4),76.4(1.0)
76-Os,,,10.7(2),,,,14.4(5),0.3(8),14.7(6),16.0(4.0)
76-Os-184,0.02,0,10.0(2.0)*,,,,13.0(5.0),0,13.0(5.0),3000.0(150.0)
76-Os-186,1.6,0,12.0(1.7),,,,17.0(5.0),0,17.0(5.0),80.0(13.0)
76-Os-187,1.6,1/2,10.0(2.0)*,,,,13.0(5.0),0.3(3)*,13.0(5.0),320.0(10.0)
76-Os-188,13.3,0,7.8(3),,,,7.3(6),0,7.3(6),4.7(5)
76-Os-189,16.1,3/2,11.0(3),,,,14.4(8),0.5(5)*,14.9(9),25.0(4.0)
76-Os-190,26.4,0,11.4(3),,,,15.2(8),0,15.2(8),13.1(3)
76-Os-192,41,0,11.9(4),,,,16.6(1.2),0,16.6(1.2),2.0(1)
77-Ir,,,10.6(3),,,,14.1(8),0.0(3.0),14.0(3.0),425.0(2.0)
77-Ir-191,37.4,3/2,12.1(9),,,,,,,954.0(10.0)
77-Ir-193,62.6,3/2,9.71(18),,,,,,,111.0(5.0)
78-Pt,,,9.60(1),,,,11.58(2),0.13(11),11.71(11),10.3(3)
78-Pt-190,0.01,0,9.0(1.0),,,,10.0(2.0),0,10.0(2.0),152.0(4.0)
78-Pt-192,1.78,0,9.9(5),,,,12.3(1.2),0,12.3(1.2),10.0(2.5)
78-Pt-194,32.9,0,10.55(8),,,,14.0(2),0,14.0(2),1.44(19)
78-Pt-195,33.8,1/2,8.91(9),9.5(3),7.2(3),+/-,9.8(2),0.13(4),9.9(2),27.5(1.2)
78-Pt-196,25.3,0,9.89(8),,,,12.3(2),0,12.3(2),0.72(4)
78-Pt-198,7.2,0,7.8(1),,,,7.6(2),0,7.6(2),3.66(19)
79-Au-197,100,3/2,7.90(7),6.26(10),9.90(14),+/-,7.32(12),0.43(5),7.75(13),98.65(9)
80-Hg,,,12.595(45),,,,20.24(5),6.6(1),26.8(1),372.3(4.0)
80-Hg-196,0.15,0,30.3(1.0),,,E,115.0(8.0),0,115.0(8.0),3080.0(180.0)
80-Hg-198,10.1,0,,,,,,0,,2.0(3)
80-Hg-199,16.9,0,16.9(4),,,E,36.0(2.0),30.0(3.0),66.0(2.0),2150.0(48.0)
80-Hg-200,23.1,0,,,,,,0,,<60.0
80-Hg-201,13.2,3/2,,,,,,,,7.8(2.0)
80-Hg-202,29.7,0,11.002(43),,,,15.2108(2),0,15.2108(2),4.89(5)
80-Hg-204,6.8,0,,,,,,0,,0.43(10)
81-Tl,,,8.776(5),,,,9.678(11),0.21(15),9.89(15),3.43(6)
81-Tl-203,29.5,1/2,8.51(8),9.08(10),6.62(10),,6.14(28),0.14(4),6.28(28),11.4(2)
81-Tl-205,70.5,1/2,8.87(7),5.15(10),9.43(10),+/-,11.39(17),0.007(1),11.40(17),0.104(17)
82-Pb,,,9.4024(13),,,,11.115(7),0.0030(7),11.118(7),0.171(2)
82-Pb-204,1.4,0,10.893(78),,,,12.3(2),0,12.3(2),0.65(7)
82-Pb-206,24.1,0,9.221(78),,,,10.68(12),0,10.68(12),0.0300(8)
82-Pb-207,22.1,1/2,9.286(16),,,+/-,10.82(9),0.002(2),10.82(9),0.699(10)
82-Pb-208,52.4,0,9.494(30),,,,11.34(5),0,11.34(5),0.00048(3)
83-Bi-209,100,9/2,8.5242(18),8.26(1),8.74(1),,9.148(4),0.0084(19),9.156(4),0.0338(7)
88-Ra-226,1620 Y,0,10.0(1.0),,,,13.0(3.0),0,13.0(3.0),12.8(1.5)
90-Th-232,100,0,10.31(3),,,,13.36(8),0,13.36(8),7.37(6)
91-Pa-231,32500 Y,3/2,9.1(3),,,,10.4(7),0.1(3.3),10.5(3.2),200.6(2.3)
92-U,,,8.417(5),,,,8.903(11),0.005(16),8.908(11),7.57(2)
92-U-233,159000 Y,5/2,10.1(2),,,,12.8(5),0.1(6),12.9(3),574.7(1.0)
92-U-234,0.005,0,12.4(3),,,,19.3(9),0,19.3(9),100.1(1.3)
92-U-235,0.72,7/2,10.50(3),,,,13.78(11),0.2(2),14.0(2),680.9(1.1)
92-U-238,99.27,0,8.407(7),,,,8.871(11),0,8.871(11),2.68(2)
93-Np-237,2140000 Y,5/2,10.55(10),,,,14.0(3),0.5(5)*,14.5(6),175.9(2.9)
94-Pu-239,24400 Y,1/2,7.7(1),,,,7.5(2),0.2(6),7.7(6),1017.3(2.1)
94-Pu-240,6540 Y,0,3.5(1),,,,1.54(9),0,1.54(9),289.6(1.4)
94-Pu-242,376000 Y,0,8.1(1),,,,8.2(2),0,8.2(2),18.5(5)
95-Am-243,7370 Y,5/2,8.3(2),,,,8.7(4),0.3(2.6),9.0(2.6),75.3(1.8)
96-Cm-244,17.9 Y,0,9.5(3),,,,11.3(7),0,11.3(7),16.2(1.2)
96-Cm-246,4700 Y,0,9.3(2),,,,10.9(5),0,10.9(5),1.36(17)
96-Cm-248,340000 Y,0,7.7(2),,,,7.5(4),0,7.5(4),3.00(26)\
"""

# Imaginary values for select isotopes
# isotope, b_c_i, bp_i, bm_i
nsftableI = """\
2-He-3,-1.48,,-5.925
3-Li-6,-0.26,-0.08(1),-0.62(2)
5-B,-0.21,,
47-Ag-107,-0.01,,
47-Ag-109,-0.025,,
48-Cd,-1.2,,
48-Cd-113,-12,,
49-In,-0.054,,
49-In-115,-0.056,,
52-Te-123,-0.1,,
62-Sm,-1.5,,
62-Sm-149,-11,,
64-Gd,-13.6,,
64-Gd-155,-10.3,,
71-Lu-176,-0.57(2),,
80-Hg-196,-0.8,,\
"""
# Excluding the following because the measurements for the real parts
# were not used in nsftable table.
# 63-Eu-151,-2.46,,
# 64-Gd-157,-47,-75,

def fix_number(str):
    """
    Converts strings of the form e.g., 35.24(2)* into numbers without
    uncertainty. Also accepts a limited range, e.g., <1e-6, which is
    converted as 1e-6.  Missing values are set to 0.
    """
    from .util import parse_uncertainty
    return parse_uncertainty(str.replace('<','').replace('*',''))[0]

def sld_table(wavelength=1, table=None, isotopes=True):
    r"""
    Scattering length density table for wavelength 4.75 |Ang|.

    :Parameters:

        *table* \: PeriodicTable
            If *table* is not specified, use the common periodic table.

        *isotopes* = True \: boolean
            Whether to consider isotopes or not.

    :Returns: None

    Example

        >>> sld_table(wavelength=4.75)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
         Neutron scattering length density table
        atom       mass density     sld    imag   incoh
        H         1.008   0.071  -1.582   0.000  10.690
        1-H       1.008   0.071  -1.582   0.000  10.691
        D         2.014   0.141   2.820   0.000   1.709
        T         3.016   0.212   2.027   0.000   0.453
        He        4.003   0.122   0.569   0.000   0.003
        3-He      3.016   0.092   1.054   0.272   0.652 *
        4-He      4.003   0.122   0.569   0.000   0.000
           ...
        248-Cm  248.072  13.569   2.536   0.000   0.207
        * Energy dependent cross sections
    """
    table = default_table(table)
    # Table for comparison with scattering length density calculators
    # b_c for Sc, Te, Xe, Sm, Eu, Gd, W, Au, Hg are different from Neutron News
    # The Rauch data have cited references to back up the numbers
    # (see doc directory), though it is not clear what criteria are
    # used to select amongst the available measurements.
    print(" Neutron scattering length density table")
    print("%-7s %7s %7s %7s %7s %7s"
          %('atom', 'mass', 'density', 'sld', 'imag', 'incoh'))
    for el in table:
        if el.neutron.has_sld():
            coh, jcoh, inc = el.neutron.sld(wavelength=wavelength)
            print("%-7s %7.3f %7.3f %7.3f %7.3f %7.3f%s"
                  %(el, el.mass, el.density, coh, jcoh, inc,
                    ' *' if el.neutron.is_energy_dependent else ''))
            if isotopes:
                isos = [iso for iso in el if iso.neutron is not None and iso.neutron.has_sld()]
            else:
                isos = []
            for iso in isos:
                coh, jcoh, inc = iso.neutron.sld(wavelength=wavelength)
                print("%-7s %7.3f %7.3f %7.3f %7.3f %7.3f%s"
                      %(iso, iso.mass, iso.density, coh, jcoh, inc,
                        ' *' if iso.neutron.is_energy_dependent else ''))
    print("* Energy dependent cross sections")

def energy_dependent_table(table=None):
    r"""
    Prints a table of energy dependent isotopes.

    :Parameters:
        *table* \: PeriodicTable
            If *table* is not specified, use the common periodic table.

    :Returns: None

    Example

        >>> energy_dependent_table()
        Elements and isotopes with energy dependent absorption:
            He-3
            Cd Cd-113
            Sm Sm-149
            Eu Eu-151
            Gd Gd-155 Gd-157
            Yb-168
            Hg-196 Hg-199
    """
    table = default_table(table)
    # List of energy dependent elements and isotopes
    print("Elements and isotopes with energy dependent absorption:")
    for el in table:
        if not hasattr(el, 'neutron'):
            continue
        dep = []
        if el.neutron.is_energy_dependent:
            dep += [str(el)]
        dep += [str(el)+'-'+str(iso.isotope)
                for iso in el
                if iso.neutron is not None and iso.neutron.is_energy_dependent]
        if dep:
            print("    " + " ".join(dep))

def _diff(iso, a, b, tol=0.01):
    if None in (a, b):
        if a is not None or b is not None:
            if a is None and b > tol:
                print("%10s %8s %8.2f"%(iso, "----", b))
            elif b is None and a > tol:
                print("%10s %8.2f %8s"%(iso, a, "----"))
    # Tricky code: Using tolerance of -tol selects for items within tolerance
    # rather than outside tolerance by using -|a-b| > -tol.
    elif np.sign(tol)*abs(a - b) > tol:
        print("%10s %8.2f %8.2f %5.1f%%"
              % (iso, a, b, (100*(a-b)/b if b != 0 else inf)))

def compare(fn1, fn2, table=None, tol=0.01):
    table = default_table(table)
    for el in table:
        try:
            res1 = fn1(el)
        except Exception:
            res1 = None
        try:
            res2 = fn2(el)
        except Exception:
            res2 = None
        _diff(el, res1, res2, tol=tol)
        for iso in el:
            # Don't show isotope details if the isotope defers to the natural
            # natural abundance for its value.
            if 'neutron' not in iso.__dict__:
                #print("dict has", iso.__dict__.keys())
                continue
            try:
                res1 = fn1(iso)
            except Exception:
                res1 = None
            try:
                res2 = fn2(iso)
            except Exception:
                res2 = None
            _diff(iso, res1, res2, tol=tol)

def absorption_comparison_table(table=None, tol=None):
    r"""
    Prints a table comparing absorption to the imaginary bound coherent
    scattering length b_c_i.  This is used to checking the integrity
    of the data and formula.

    The relationship between absorption and b_c_i is:

    .. math::

        \sigma_a = -2 \lambda \mathrm{Im}(b_c) \cdot 1000

    The wavelength $\lambda = 1.798$ |Ang| is the neutron wavelength at which
    the absorption is tallied. The factor of 1000 transforms from
    |Ang|\ |cdot|\ fm to barn.

    :Parameters:
        *table* \: PeriodicTable
            The default periodictable unless a specific table has been requested.
        *tol* = 0.01 \: float | barn
            Show differences greater than this amount.

    :Returns: None

    Example

        >>> absorption_comparison_table (tol=0.5) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        Comparison of absorption and (-2000 lambda b_c_i)
              3-He  5333.00  5322.08   0.2%
                Li    70.50     ----
              6-Li   940.00   934.96   0.5%
                 B   767.00   755.16   1.6%
              10-B  3835.00     ----
                 N     1.90     ----
           ...

    """

    print("Comparison of absorption and (-2000 lambda b_c_i)")
    compare(lambda el: el.neutron.absorption,
            lambda el: -2000*el.neutron.b_c_i*ABSORPTION_WAVELENGTH,
            table=table, tol=tol)

def coherent_comparison_table(table=None, tol=None):
    r"""
    Prints a table of $4 \pi |b_c|^2/100$ and coherent for each isotope.
    This is useful for checking the integrity of the data and formula.

    The table only prints where b_c exists.

    :Parameters:
        *table* \: PeriodicTable
            The default periodictable unless a specific table has been requested.
        *tol* = 0.01 \: float | barn
            Amount of difference to show. Use -tol to show elements within
            tolerance rather than those outside tolerance.

    :Returns: None

    Example

        >>> coherent_comparison_table (tol=0.5) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        Comparison of (4 pi |b_c|^2/100) and coherent
                Sc    18.40    19.00  -3.2%
             45-Sc    18.40    19.00  -3.2%
             65-Cu    13.08    14.10  -7.2%
             84-Sr     3.14     6.00 -47.6%
           ...

    """
    print("Comparison of (4 pi |b_c|^2/100) and coherent")
    sigma_c = lambda el: 4*pi/100*abs(el.neutron.b_c_complex)**2
    compare(sigma_c, lambda el: el.neutron.coherent, table=table, tol=tol)

def total_comparison_table(table=None, tol=None):
    r"""
    Prints a table of neutron.total and sum coh,inc for each
    isotope where these exist.  This is used to checking the integrity
    of the data and formula.

    :Parameters:
        *table* \: PeriodicTable
            The default periodictable unless a specific table has been requested.
        *tol* = 0.01 \: float | barn
            Amount of difference to show. Use -tol to show elements within
            tolerance rather than those outside tolerance.

    :Returns: None

    Example

        >>> total_comparison_table (tol=0.1)
        Comparison of total cross section to (coherent + incoherent)
             84-Kr     6.60     ----
                Xe     4.34     3.04  42.9%
            149-Sm   200.00   200.50  -0.2%
                Eu     9.20     9.07   1.4%
                Gd   180.00   180.30  -0.2%
            155-Gd    66.00    65.80   0.3%
            161-Dy    16.00    16.30  -1.8%
            180-Ta     7.00     6.70   4.5%
            187-Os    13.00    13.30  -2.3%

    """
    print("Comparison of total cross section to (coherent + incoherent)")
    compare(lambda el: el.neutron.total,
            lambda el: el.neutron.coherent+el.neutron.incoherent,
            table=table, tol=tol)

def incoherent_comparison_table(table=None, tol=None):
    r"""
    Prints a table of incoherent computed from total and b_c with incoherent.

    :Parameters:
        *table* \: PeriodicTable
            The default periodictable unless a specific table has been requested.
        *tol* = 0.01 \: float | barn
            Amount of difference to show. Use -tol to show elements within
            tolerance rather than those outside tolerance.

    :Returns: None

    Example

        >>> incoherent_comparison_table (tol=0.5) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        Comparison of incoherent and (total - 4 pi |b_c|^2/100)
                Sc     4.50     5.10 -11.8%
             45-Sc     4.50     5.10 -11.8%
             65-Cu     0.40     1.42 -71.7%
             84-Sr     0.00     2.86 -100.0%
           ...

    """
    print("Comparison of incoherent and (total - 4 pi |b_c|^2/100)")
    sigma_c = lambda el: 4*pi/100*abs(el.neutron.b_c_complex)**2
    compare(lambda el: el.neutron.incoherent,
            lambda el: el.neutron.total - sigma_c(el),
            table=table, tol=tol)

def print_scattering(compound, wavelength=ABSORPTION_WAVELENGTH):
    """
    Print the scattering for a single compound.
    """
    from . import formulas
    compound = formulas.formula(compound)
    density = compound.density if compound.density is not None else 1.0
    sld, xs, penetration = neutron_scattering(compound, wavelength=wavelength,
                                              density=density)
    print("%s at %g Ang  (density=%g g/cm^3)"
          % (str(compound), wavelength, density))
    print("  sld: %g + %g j  (%g incoherent)  1e-6/Ang^2"%sld)
    print("  Σ_c: %g  Σ_a: %g  Σ_i: %g  1/cm"%xs)
    print("  μ: %g 1/cm  1/e penetration: %g cm"%(1/penetration, penetration))

def main():
    """
    Simple command line interface, showing the predicted neutron scattering.

    Usage::

        python -m periodictable.nsf [-Lwavelength] compound@density compound@density ...

    For example::

        $ python -m periodictable.nsf XeF6@3.56
        scattering for XeF6 at 1.798 Ang  (density=3.56 g/cm^3)
          sld: 3.37503 + 0.000582313 j  (0.402605 incoherent)  1e-6/Ang^2
          sigma_c: 3.37503  sigma_i: 0.000582313  sigma_a: 0.402605  1/cm
          1/e penetration: 2.23871 cm
    """

    import sys
    compounds = sys.argv[1:]
    if compounds[0].startswith('-L'):
        wavelength = float(compounds[0][2:])
        compounds = compounds[1:]
    else:
        wavelength = ABSORPTION_WAVELENGTH
    for c in compounds:
        print_scattering(c, wavelength)

if __name__ == "__main__":
    main()
    #sld_table()
    #coherent_comparison_table(tol=0.1)
    #incoherent_comparison_table(tol=0.1)
    #absorption_comparison_table(tol=0.1)
    #total_comparison_table(tol=0.1)
