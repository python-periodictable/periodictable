.. _using:

***********
Basic usage
***********

The periodic table is available on `PyPI <http://pypi.python.org/pypi>`_,
and can be obtained simply with::

    easy_install periodictable

This will install pyparsing if it is not already available.  The numpy
package must already be installed.

Access particular elements by name:

.. doctest::

    >>> from periodictable import hydrogen
    >>> print("H mass %s %s"%(hydrogen.mass, hydrogen.mass_units))
    H mass 1.008 u

Access particular elements as symbols:

.. doctest::

    >>> from periodictable import H,B,Cu,Ni
    >>> print("B absorption %s"%B.neutron.absorption)
    B absorption 767.0
    >>> print("Ni f1/f2 for Cu K-alpha X-rays f'=%.5f f''=%.5f"
    ...       % Ni.xray.scattering_factors(wavelength=Cu.K_alpha))
    Ni f1/f2 for Cu K-alpha X-rays f'=25.02293 f''=0.52491

Access isotopes using mass number subscripts:

.. doctest::

    >>> print("58-Ni vs 62-Ni scattering %s:%s"%(Ni[58].neutron.coherent, Ni[62].neutron.coherent))
    58-Ni vs 62-Ni scattering 26.1:9.5

Access elements indirectly:

.. doctest::

    >>> import periodictable
    >>> print("Cd density %.2f %s"%(periodictable.Cd.density, periodictable.Cd.density_units))
    Cd density 8.65 g/cm^3

Import all elements:

.. doctest::

    >>> from periodictable import *
    >>> print(periodictable.H)
    H
    >>> print(periodictable.H.mass)
    1.008

Deuterium and tritium are special isotopes named D and T:

.. doctest::

    >>> print("D mass %s"%D.mass)
    D mass 2.01410177784

Process all the elements:

.. doctest::

    >>> import periodictable
    >>> for el in periodictable.elements: # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    ...     print("%s %s"%(el.symbol,el.name))
    H hydrogen
    He helium
       ...
    Og oganesson

Another example for processing all elements:

.. doctest::

    >>> from periodictable import elements
    >>> for el in elements: # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    ...     print("%s %s"%(el.symbol,el.number))
    H 1
    He 2
       ...

Process all the :class:`isotopes <periodictable.core.Isotope>` for an element:

.. doctest::

    >>> for iso in periodictable.H:
    ...     print("%s %s"%(iso,iso.mass))
    1-H 1.0078250319
    D 2.01410177784
    T 3.01604928132
    4-H 4.02643
    5-H 5.03531
    6-H 6.04496

You can create a unique handle to an individual ion.  In addition to storing
the ion charge, this can be used to reference the underlying properties of
the element or isotope:

.. doctest::

    >>> Ni58_2 = periodictable.Ni[58].ion[2]
    >>> Ni_2 = periodictable.Ni.ion[2]
    >>> print("charge for Ni2+ is %d"%Ni_2.charge)
    charge for Ni2+ is 2
    >>> print("mass for Ni[58] and for natural abundance: %.4f %.4f"%(Ni58_2.mass, Ni_2.mass))
    mass for Ni[58] and for natural abundance: 57.9342 58.6923

The ion specific properties can be accessed from the ion using ion.charge
for the ion index:

.. doctest::

    >>> import periodictable
    >>> Fe_2 = periodictable.Fe.ion[2]
    >>> print("[%.5f, %.5f, %.5f]"
    ...       % tuple(Fe_2.magnetic_ff[Fe_2.charge].M_Q([0,0.1,0.2])))
    [1.00000, 0.99935, 0.99741]

The following is a plot of the magnetic form factor vs. Q:

.. doctest::

    >>> import pylab # doctest: +SKIP
    >>> Q = pylab.linspace(0,16,200) # doctest: +SKIP
    >>> M = Fe_2.magnetic_ff[Fe_2.charge].j0_Q(Q) # doctest: +SKIP
    >>> pylab.xlabel(r'Magnetic Form Factor for Fe') # doctest: +SKIP
    >>> pylab.ylabel(r'$\AA^{-1}$') # doctest: +SKIP
    >>> pylab.title('Ion specific property for Fe') # doctest: +SKIP
    >>> pylab.plot(Q,M) # doctest: +SKIP

.. plot:: plots/magnetic_ff.py

Missing properties generally evaluate to *None*:

.. doctest::

    >>> print("Radon density %s"%periodictable.Rn.density)
    Radon density None


Specific defined properties related to elements can be accessed in a table format as shown in following example :

.. doctest::

    >>> elements.list('symbol','K_alpha',format="%s K-alpha = %s") # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Ne K-alpha = 14.6102
    Na K-alpha = 11.9103
    Mg K-alpha = 9.8902
    Al K-alpha = 8.3402
       ...
    Cf K-alpha = 0.1094
    Es K-alpha = 0.1067
    Fm K-alpha = 0.104
