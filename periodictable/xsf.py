# This program is public domain
# Author: Paul Kienzle
r"""
This module has one class and nine fuctions.

   :class:`Xray`
       X-ray scattering properties for the elements.

The following attributes are added to each element:

   :func:`Xray.sftable`
       Three column table of energy vs. scattering factors f1, f2.

   :func:`Xray.scattering_factors`
       Returns f1, f2, the X-ray scattering factors for the given wavelengths
       interpolated from sftable.

   :func:`Xray.f0`
       Returns f0 for the given vector Q, with Q[i] in $[0, 24\pi]$ |1/Ang|.

   :func:`Xray.sld`
       Returns scattering length density (*real*, *imaginary*) for the
       given wavelengths or energies.

The following functions are available for X-ray scatting information processing:

     :func:`xray_wavelength`
         Finds X-ray wavelength in angstroms given energy in keV.

     :func:`xray_energy`
         Finds X-ray energy in keV given wavelength in angstroms.

     :func:`init`
         Initializes a periodic table with the Lawrence Berkeley Laboratory
         Center for X-Ray Optics xray scattering factors.

     :func:`init_spectral_lines`
         Sets the K_alpha and K_beta1 wavelengths for select elements.

     :func:`sld_table`
         Prints the xray SLD table for the given wavelength.

     :func:`xray_sld`
         Computes xray scattering length densities for molecules.

     :func:`index_of_refraction`
         Express xray scattering length density as an index of refraction

     :func:`mirror_reflectivity`
         X-ray reflectivity from a mirror made of a single compound.

     :func:`xray_sld_from_atoms`
         The underlying scattering length density calculator. This works with
         a dictionary of atoms and quantities directly.

     :func:`emission_table`
         Prints a table of emission lines.

K_alpha, K_beta1 (|Ang|):
    X-ray emission lines for elements beyond neon, with
    $K_\alpha = (2 K_{\alpha 1} + K_{\alpha 2})/3$.

X-ray scattering factors:
    Low-Energy X-ray Interaction Coefficients: Photoabsorption, scattering
    and reflection for E in 30 to 30,000 eV, and Z in 1 to 92.

.. Note::

    For :ref:`custom tables <custom-table>`, use :func:`init` and
    :func:`init_spectral_lines` to set the data.


Emission line tables
====================

Data for the $K_\alpha$ and $K_\beta$ lines comes from
[#Deslattes2003]_, with the full tables available at
`<http://www.nist.gov/pml/data/xraytrans/index.cfm>`_.
Experimental Values are used, truncated to 4 digits
of precision to correspond to the values for the subset
of elements previously defined in the periodictable package.

X-ray f1 and f2 tables
======================
The data for the tables is stored in the ``periodictable/xsf``.
directory.  The following information is from ``periodictable/xsf/read.me``,
with minor formatting changes:

  These ``[*.nff]`` files were used to generate the tables published in
  reference [#Henke1993]_. The files contain three columns of data:

    Energy(eV), *f_1*, *f_2*,

  where *f_1* and *f_2* are the atomic (forward) scattering factors.
  There are 500+ points on a uniform logarithmic mesh with points
  added 0.1 eV above and below "sharp" absorption edges. The
  tabulated values of *f_1* contain a relativistic, energy
  independent, correction given by:

  .. math::

    Z^* = Z - (Z/82.5)^{2.37}

  .. Note::
    Below 29 eV *f_1* is set equal to -9999.

  The atomic photoabsorption cross section, $\mu_a$, may be readily
  obtained from the values of $f_2$ using the relation:

  .. math::

    \mu_a = 2 r_e \lambda f_2

  where $r_e$ is the classical electron radius, and $\lambda$ is
  the wavelength. The index of refraction for a material with *N* atoms per
  unit volume is calculated by:

  .. math::

    n = 1 - N r_e \lambda^2 (f_1 + i f_2)/(2 \pi).

  These (semi-empirical) atomic scattering factors are based upon
  photoabsorption measurements of elements in their elemental state.
  The basic assumption is that condensed matter may be modeled as a
  collection of non-interacting atoms.  This assumption is in general
  a good one for energies sufficiently far from absorption thresholds.
  In the threshold regions, the specific chemical state is important
  and direct experimental measurements must be made.

  These tables are based on a compilation of the available experimental
  measurements and theoretical calculations.  For many elements there is
  little or no published data and in such cases it was necessary to
  rely on theoretical calculations and interpolations across Z.
  In order to improve the accuracy in the future considerably more
  experimental measurements are needed.

  Note that the following elements have been updated since the
  publication of Ref. [#Henke1993]_ in July 1993. Check
  `<http://henke.lbl.gov/optical_constants/update.html>`_ for more
  recent updates.

  .. table::

    ========  ==========   =================
    Element   Updated      Energy Range (eV)
    ========  ==========   =================
    Pt        Aug 2023     0.1-800
    Cr        Aug 2023     0.4-813
    Nb        Aug 2023     10-1000
    Y         Aug 2023     50-1300
    Er        Aug 2023     3-1580
    Mg        Jan 2011     10-1300
    Zr        Apr 2010     20-1000
    La        Jun 2007     14-440
    Gd        Jun 2007     12-450
    Sc        Apr 2006     50-1300
    Ti        Aug 2004     20-150
    Ru        Aug 2004     40-1300
    W         Aug 2004     35-250
    Mo        Aug 2004     25-60
    Be        Aug 2004     40-250
    Mo        Nov 1997     10-930
    Fe        Oct 1995     600-800
    Si        Jun 1995     30-500
    Au        Jul 1994     2000-6500
    Mg,Al,Si  Jan 1994     30-200
    Li        Nov 1994     2000-30000
    ========  ==========   =================


  Data available at:
    #. http://henke.lbl.gov/optical_constants/asf.html

.. [#Henke1993] B. L. Henke, E. M. Gullikson, and J. C. Davis.  "X-ray interactions:
       photoabsorption, scattering, transmission, and reflection at E=50-30000 eV,
       Z=1-92", Atomic Data and Nuclear Data Tables 54 no.2, 181-342 (July 1993).

.. [#Deslattes2003] R. D. Deslattes, E. G. Kessler, Jr., P. Indelicato, L. de Billy,
       E. Lindroth, and J. Anton.  Rev. Mod. Phys. 75, 35-99 (2003).
"""
__all__ = ['Xray', 'init', 'init_spectral_lines',
           'xray_energy', 'xray_wavelength',
           'xray_sld', 'xray_sld_from_atoms',
           'emission_table', 'sld_table', 'plot_xsf',
           'index_of_refraction', 'mirror_reflectivity',
           ]
import os.path

import numpy as np
from numpy import nan, pi, exp, sin, cos, sqrt, radians

from .core import Element, Ion, default_table, get_data_path
from .constants import (
    avogadro_number, planck_constant, speed_of_light, electron_volt,
    electron_radius)

def xray_wavelength(energy):
    r"""
    Convert X-ray energy to wavelength.

    :Parameters:
        *energy* : float or vector | keV

    :Returns:
        *wavelength* : float | |Ang|

    Energy can be converted to wavelength using

    .. math::

        \lambda = h c / E

    where:

        $h$ = Planck constant in J\ |cdot|\ s

        $c$ = speed of light in m/s
    """
    return planck_constant/electron_volt*speed_of_light/np.asarray(energy)*1e7

def xray_energy(wavelength):
    r"""
    Convert X-ray wavelength to energy.

    :Parameters:
        *wavelength* : float or vector | |Ang|

    :Returns:
        *energy* : float or vector | keV

    Wavelength can be converted to energy using

    .. math::

        E = h c / \lambda

    where:

        $h$ = Planck constant in J\ |cdot|\ s

        $c$ = speed of light in m/s
    """
    return planck_constant/electron_volt*speed_of_light/np.asarray(wavelength)*1e7

class Xray:
    """
    X-ray scattering properties for the elements. Refer help(periodictable.xsf)
    from command prompt for details.
    """
    sftable_units: tuple[str, str, str] = ("eV", "", "")
    scattering_factors_units: tuple[str, str] = ("", "")
    sld_units: tuple[str, str] = ("1e-6/Ang^2", "1e-6/Ang^2")
    element: Element
    _nff_path = get_data_path('xsf')
    _table = None
    def __init__(self, element):
        self.element = element

    def _gettable(self):
        if self._table is None:
            # Load table when necessary; note there is no table for
            # neutrons (n), and lowercase nitrogen=> n.nff, so it must
            # be checked for explicitly.
            filename = os.path.join(self._nff_path,
                                    self.element.symbol.lower()+".nff")
            if self.element.symbol != 'n' and os.path.exists(filename):
                xsf = np.loadtxt(filename, skiprows=1).T
                xsf[1, xsf[1] == -9999.] = np.nan
                xsf[0] *= 0.001  # Use keV in table rather than eV
                self._table = xsf
        return self._table
    sftable = property(_gettable, doc="X-ray scattering factor table (E,f1,f2)")

    def scattering_factors(self, *, energy=None, wavelength=None):
        """
        X-ray scattering factors f', f''.

        :Parameters:
            *energy* : float or vector | keV
                X-ray energy.

        :Returns:
            *scattering_factors* : (float, float)
                Values outside the range return NaN.

        Values are found from linear interpolation within the Henke Xray
        scattering factors database at the Lawrence Berkeley Laboratory
        Center for X-ray Optics.
        """
        xsf = self.sftable
        if xsf is None:
            return None, None

        if wavelength is not None:
            energy = xray_energy(wavelength)
        if energy is None:
            raise TypeError('X-ray scattering factors need wavelength or energy')

        scalar = np.isscalar(energy)
        if scalar:
            energy = np.array([energy])
        f1 = np.interp(energy, xsf[0], xsf[1], left=nan, right=nan)
        f2 = np.exp(np.interp(
            np.log(energy), np.log(xsf[0]), np.log(xsf[2]),
            left=nan, right=nan))
        if scalar:
            f1, f2 = f1[0], f2[0]
        return f1, f2

    def f0(self, Q):
        r"""
        Isotropic X-ray scattering factors *f0* for the input Q.

        :Parameters:
            *Q* : float or vector in $[0, 24\pi]$ | |1/Ang|
                X-ray scattering properties for the elements.

        :Returns:
            *f0* : float
                Values outside the valid range return NaN.


        .. Note::

            *f0* is often given as a function of $\sin(\theta)/\lambda$
            whereas we are using  $Q = 4 \pi \sin(\theta)/\lambda$, or
            in terms of energy $Q = 4 \pi \sin(\theta) E/(h c)$.

        Reference:
             D. Wassmaier, A. Kerfel, Acta Crystallogr. A51 (1995) 416.
             http://dx.doi.org/10.1107/S0108767394013292
        """
        from . import cromermann
        f = cromermann.fxrayatq(Q=Q,
                                symbol=self.element.symbol,
                                charge=self.element.charge)
        return f

    def sld(self, *, wavelength=None, energy=None):
        r"""
        X ray scattering length density.

        :Parameters:
            *wavelength* : float or vector | |Ang|
                Wavelength of the X-ray.

            *energy* : float or vector | keV
                Energy of the X-ray (if *wavelength* not specified).

            .. note:
                Only one of *wavelength* and *energy* is needed.

        :Returns:
            *sld* : (float, float) | |1/Ang^2|
                (*real*, *imaginary*) X-ray scattering length density.

        :Raises:
            *TypeError* : neither *wavelength* nor *energy* was specified.

        The scattering length density is $r_e N (f_1 + i f_2)$.
        where $r_e$ is the electron radius and $N$ is the
        number density.  The number density is $N = \rho_m/m N_A$,
        with mass density $\rho_m$ molar mass $m$ and
        Avogadro's number $N_A$.

        The constants are available directly:

            $r_e$ = periodictable.xsf.electron_radius

            $N_A$ = periodictable.constants.avogadro_number

        Data comes from the Henke Xray scattering factors database at the
        Lawrence Berkeley Laboratory Center for X-ray Optics.
        """
        f1, f2 = self.scattering_factors(wavelength=wavelength, energy=energy)
        if f1 is None or self.element.number_density is None:
            return None, None
        rho = f1*electron_radius*self.element.number_density*1e-8
        irho = f2*electron_radius*self.element.number_density*1e-8
        return rho, irho

# Note: docs and function prototype are reproduced in __init__
def xray_sld(compound, *, density=None, natural_density=None,
             wavelength=None, energy=None):
    """
    Compute xray scattering length densities for molecules.

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.
        *density* : float | |g/cm^3|
            Mass density of the compound, or None for default.
        *natural_density* : float | |g/cm^3|
            Mass density of the compound at naturally occurring isotope abundance.
        *wavelength* : float or vector | |Ang|
            Wavelength of the X-ray.
        *energy* : float or vector | keV
            Energy of the X-ray, if *wavelength* is not specified.

    :Returns:
        *sld* : (float, float) | |1e-6/Ang^2|
            (*real*, *imaginary*) scattering length density.

    :Raises:
        *AssertionError* :  *density* or *wavelength*/*energy* is missing.
    """
    from . import formulas
    compound = formulas.formula(compound, density=density,
                                natural_density=natural_density)
    assert compound.density is not None, "scattering calculation needs density"

    if wavelength is not None:
        energy = xray_energy(wavelength)
    assert energy is not None, "scattering calculation needs energy or wavelength"

    mass, sum_f1, sum_f2 = 0, 0, 0
    for element, quantity in compound.atoms.items():
        mass += element.mass*quantity
        f1, f2 = element.xray.scattering_factors(energy=energy)
        if f1 is None:
            raise ValueError(f"X-ray scattering factors not available for {element}")
        #print element, f1, f2, wavelength
        sum_f1 += f1*quantity
        sum_f2 += f2*quantity

    if mass == 0: # because the formula is empty
        return 0, 0

    N = (compound.density/mass*avogadro_number*1e-8)
    rho = N*sum_f1*electron_radius
    irho = N*sum_f2*electron_radius
    return rho, irho


def index_of_refraction(compound, *, density=None, natural_density=None,
                        energy=None, wavelength=None):
    """
    Calculates the index of refraction for a given compound

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.
        *density* : float | |g/cm^3|
            Mass density of the compound, or None for default.
        *natural_density* : float | |g/cm^3|
            Mass density of the compound at naturally occurring isotope abundance.
        *wavelength* : float or vector | |Ang|
            Wavelength of the X-ray.
        *energy* : float or vector | keV
            Energy of the X-ray, if *wavelength* is not specified.

    :Returns:
        *n* : float or vector | unitless
            index of refraction of the material at the given energy

    :Notes:

    Formula taken from http://xdb.lbl.gov (section 1.7) and checked
    against http://henke.lbl.gov/optical_constants/getdb2.html
    """
    if energy is not None:
        wavelength = xray_wavelength(energy)
    assert wavelength is not None, "scattering calculation needs energy or wavelength"
    f1, f2 = xray_sld(compound,
                      density=density, natural_density=natural_density,
                      wavelength=wavelength)
    return 1 - wavelength**2/(2*pi)*(f1 + f2*1j)*1e-6

def delta(compound, *, density=None, natural_density=None,
                        energy=None, wavelength=None):
    """
    Calculates the δ component of the index of refraction for a given compound

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.
        *density* : float | |g/cm^3|
            Mass density of the compound, or None for default.
        *natural_density* : float | |g/cm^3|
            Mass density of the compound at naturally occurring isotope abundance.
        *wavelength* : float or vector | |Ang|
            Wavelength of the X-ray.
        *energy* : float or vector | keV
            Energy of the X-ray, if *wavelength* is not specified.

    :Returns:
        *delta* : float or vector | unitless
            δ component of the index of refraction of the material at the given energy

    :Notes:

    Formula taken from http://xdb.lbl.gov (section 1.7) and checked
    against http://henke.lbl.gov/optical_constants/getdb2.html.
    """
    if energy is not None:
        wavelength = xray_wavelength(energy)
    assert wavelength is not None, "scattering calculation needs energy or wavelength"
    f1, f2 = xray_sld(compound,
                      density=density, natural_density=natural_density,
                      wavelength=wavelength)
    return wavelength**2/(2*pi)*(f1)*1e-6

def beta(compound, *, density=None, natural_density=None,
                        energy=None, wavelength=None):
    r"""
    Calculates the β component of the index of refraction for a given compound

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.
        *density* : float | |g/cm^3|
            Mass density of the compound, or None for default.
        *natural_density* : float | |g/cm^3|
            Mass density of the compound at naturally occurring isotope abundance.
        *wavelength* : float or vector | |Ang|
            Wavelength of the X-ray.
        *energy* : float or vector | keV
            Energy of the X-ray, if *wavelength* is not specified.

    :Returns:
        *beta* : float or vector | unitless
            β component of the index of refraction of the material at the given energy

    :Notes:

    Formula taken from the equation

    .. math::

        n = 1 - \delta + i \beta

    derivative of the formula http://xdb.lbl.gov (section 1.7) using a more common
    sign convention for the imaginary part of the index of refraction. Checked
    against http://henke.lbl.gov/optical_constants/getdb2.html.
    """
    if energy is not None:
        wavelength = xray_wavelength(energy)
    assert wavelength is not None, "scattering calculation needs energy or wavelength"
    f1, f2 = xray_sld(compound,
                      density=density, natural_density=natural_density,
                      wavelength=wavelength)
    return wavelength**2/(2*pi)*(f2)*1e-6

def mirror_reflectivity(compound, *, density=None, natural_density=None,
                        energy=None, wavelength=None,
                        angle=None, roughness=0):
    """
    Calculates reflectivity of a thick mirror as function of energy and angle

    :Parameters:
        *compound* : Formula initializer
            Chemical formula.
        *density* : float | |g/cm^3|
            Mass density of the compound, or None for default.
        *natural_density* : float | |g/cm^3|
            Mass density of the compound at naturally occurring isotope abundance.
        *roughness* : float | |Ang|
            High-spatial-frequency surface roughness.
        *wavelength* : float or vector | |Ang|
            Wavelength of the X-ray.
        *energy* : float or vector | keV
            Energy of the X-ray, if *wavelength* is not specified.
        *angle* : vector | |deg|
            Incident beam angles.

    :Returns:
        *reflectivity* : matrix
            matrix of reflectivity as function of (angle, energy)

    :Notes:

    Formula taken from http://xdb.lbl.gov (section 4.2) and checked
    against http://henke.lbl.gov/optical_constants/mirror2.html
    """
    if energy is not None:
        wavelength = xray_wavelength(energy)
    assert wavelength is not None, "scattering calculation needs energy or wavelength"
    assert angle is not None, "scattering calculation needs incident angle"
    angle = radians(angle)
    if np.isscalar(wavelength):
        wavelength = np.array([wavelength])
    if np.isscalar(angle):
        angle = np.array([angle])
    nv = index_of_refraction(compound=compound,
                             density=density, natural_density=natural_density,
                             wavelength=wavelength)
    ki = 2*pi/wavelength[None, :] * sin(angle[:, None])
    kf = 2*pi/wavelength[None, :] * sqrt(nv[None, :]**2 - cos(angle[:, None])**2)
    r = (ki-kf)/(ki+kf)*exp(-2*ki*kf*roughness**2)
    return abs(r)**2


def xray_sld_from_atoms(*args, **kw):
    """
    .. deprecated:: 0.91

        :func:`xray_sld` now accepts a dictionary of *{atom: count}* directly.
    """
    return xray_sld(*args, **kw)


spectral_lines_data = """\
Ac  0.1380  0.1205
Ag  0.5608  0.4970
Al  8.3402  7.9601
Am  0.1181  0.1030
Ar  4.1929  3.8860
As  1.1772  1.0573
At  0.1537  0.1343
Au  0.1818  0.1589
Ba  0.3866  0.3408
Bi  0.1624  0.1419
Bk  0.1122  0.0979
Br  1.0411  0.9328
Ca  3.3595  3.0897
Cd  0.5364  0.4751
Ce  0.3586  0.3158
Cf  0.1094  0.0954
Cl  4.7287  4.4034
Cm  0.1151  0.1004
Co  1.7902  1.6208
Cr  2.2910  2.0848
Cs  0.4018  0.3543
Cu  1.5418  1.3922
Dy  0.2711  0.2378
Er  0.2539  0.2226
Es  0.1067  0.0930
Eu  0.3000  0.2635
Fe  1.9373  1.7566
Fm  0.1040  0.0907
Fr  0.1456  0.1271
Ga  1.3414  1.2079
Gd  0.2899  0.2546
Ge  1.2553  1.1289
Hf  0.2238  0.1960
Hg  0.1767  0.1544
Ho  0.2623  0.2301
I   0.4348  0.3839
In  0.5136  0.4545
Ir  0.1926  0.1685
K   3.7423  3.4539
Kr  0.9816  0.8785
La  0.3722  0.3279
Lu  0.2308  0.2023
Mg  9.8902  9.5211
Mn  2.1031  1.9102
Mo  0.7107  0.6323
Na 11.9103 11.5752
Nb  0.7476  0.6657
Nd  0.3333  0.2933
Ne 14.6102 14.4522
Ni  1.6592  1.5001
Np  0.1243  0.1085
Os  0.1984  0.1736
P   6.1581  5.7961
Pa  0.1310  0.1143
Pb  0.1670  0.1459
Pd  0.5869  0.5205
Pm  0.3217  0.2829
Po  0.1580  0.1380
Pr  0.3456  0.3042
Pt  0.1871  0.1636
Pu  0.1212  0.1057
Ra  0.1417  0.1238
Rb  0.9269  0.8286
Re  0.2043  0.1788
Rh  0.6147  0.5456
Rn  0.1496  0.1306
Ru  0.6445  0.5724
S   5.3731  5.0316
Sb  0.4718  0.4170
Sc  3.0320  2.7796
Se  1.1061  0.9921
Si  7.1263  6.7531
Sm  0.3105  0.2730
Sn  0.4920  0.4352
Sr  0.8766  0.7829
Ta  0.2171  0.1900
Tb  0.2802  0.2460
Tc  0.6764  0.6013
Te  0.4527  0.4000
Th  0.1344  0.1174
Ti  2.7497  2.5139
Tl  0.1717  0.1501
Tm  0.2459  0.2155
U   0.1276  0.1114
V   2.5048  2.2844
W   0.2106  0.1843
Xe  0.4178  0.3687
Y   0.8302  0.7407
Yb  0.2382  0.2088
Zn  1.4364  1.2952
Zr  0.7873  0.7018\
"""

def init_spectral_lines(table):
    """
    Sets the K_alpha and K_beta1 wavelengths for select elements
    """
    Element.K_alpha_units = "angstrom"
    Element.K_beta1_units = "angstrom"
    for row in spectral_lines_data.split('\n'):
        el, K_alpha, K_beta1 = row.split()
        el = table.symbol(el)
        el.K_alpha = float(K_alpha)
        el.K_beta1 = float(K_beta1)

def init(table, reload=False):

    if 'xray' in table.properties and not reload:
        return
    table.properties.append('xray')

    # Create an xray object for the particular element/ion.  Note that
    # we must not use normal attribute tests such as "hasattr(el, 'attr')"
    # or "try: el.attr; except:" since the delegation methods on Ion will
    # just return the attribute from the base element.  Instead we check
    # for an instance specific xray object for the particular ion prior
    # to delegating.
    # TODO: is there a better way to set up delegation on a field by
    # field basis?
    def _cache_xray(el):
        if '_xray' not in el.__dict__ and isinstance(el, (Element, Ion)):
            el._xray = Xray(el)
        return el._xray
    Element.xray = property(_cache_xray)
    Ion.xray = property(_cache_xray)

    ## Note: the simple approach below fails for e.g., Ni[58].ion[3].xray
    #for el in table:
    #    for charge in el.ions:
    #        el.ion[charge].xray = Xray(el.ion[charge])
    #    el.xray = Xray(el)

def plot_xsf(el):
    """
    Plots the xray scattering factors for the given element.

    :Parameters:
        *el* : Element

    :Returns: None
    """
    import matplotlib.pyplot as plt

    xsf = el.xray.sftable
    plt.title('X-ray scattering factors for '+el.name)
    plt.plot(xsf[0], xsf[1])
    plt.plot(xsf[0], xsf[2])
    plt.xlabel('Energy (keV)')
    plt.ylabel('Scattering factor')
    plt.legend(['f1', 'f2'])
    plt.show()

def sld_table(wavelength=None, table=None):
    """
    Prints the xray SLD table for the given wavelength.

    :Parameters:
        *wavelength* = Cu K-alpha : float | |Ang|
            X-ray wavelength.
        *table* : PeriodicTable
            The default periodictable unless a specific table has been requested.

    :Returns: None

    Example

        >>> sld_table() # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
        X-ray scattering length density for 1.5418 Ang
         El    rho   irho
          H   1.19   0.00
         He   1.03   0.00
         Li   3.92   0.00
         Be  13.93   0.01
          B  18.40   0.01
          C  18.71   0.03
          N   6.88   0.02
          O   9.74   0.04
          F  12.16   0.07
         Ne  10.26   0.09
         Na   7.98   0.09
         Mg  14.78   0.22
          ...
    """
    table = default_table(table)
    if wavelength is None:
        wavelength = table.Cu.K_alpha

    # NBCU spreadsheet format
    print("X-ray scattering length density for %g Ang"%wavelength)
    print("%3s %6s %6s"%('El', 'rho', 'irho'))
    for el in table:
        rho, irho = el.xray.sld(wavelength=table.Cu.K_alpha)
        if rho is not None:
            print("%3s %6.2f %6.2f"%(el.symbol, rho, irho))

def emission_table(table=None):
    """
    Prints a table of emission lines.

    :Parameters:
        *table* : PeriodicTable
            The default periodictable unless a specific table has been requested.

    :Returns: None

    Example

        >>> emission_table() # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
         El  Kalpha  Kbeta1
         Ne 14.6102 14.4522
         Na 11.9103 11.5752
         Mg  9.8902  9.5211
         Al  8.3402  7.9601
         Si  7.1263  6.7531
         ...
    """
    table = default_table(table)
    print("%3s %7s %7s"%('El', 'Kalpha', 'Kbeta1'))
    for el in table:
        if hasattr(el, 'K_alpha'):
            print("%3s %7.4f %7.4f"%(el.symbol, el.K_alpha, el.K_beta1))
