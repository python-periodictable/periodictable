=========================================
Extensible periodic table of the elements
=========================================

This package provides a periodic table of the elements with
support for mass, density and xray/neutron scattering information.

Neutron scattering calculations use values collected by the
Atomic Institute of the Austrian Universities as they appear in the neutron
data booklet, with support for some energy dependent scattering
in rare earth elements given by Lynn and Seeger (1990). X-ray scattering
calculations use a combination of empirical and theoretical values from
the LBL Center for X-ray Optics.

Tabulated values differ from those given in other sources such as the
International Tables for Crystallography, Volume C, and so computed
cross sections may give different results from other packages.

Neutron activation calculations are based on Shleien (1998), with
isotopes important to health physics. They do not perform a full
activation analysis, but instead give a gross estimate of the amount
of activation expected for a sample in the beam.

Install using::

    pip install periodictable

There are several web interfaces to this package:

* `NIST Center for Neutron Research (NCNR) <https://www.ncnr.nist.gov/resources/activation>`_
* `Forschungs-Neutronenquelle Heinz Maier-Leibnitz (FRM II) <https://webapps.frm2.tum.de/intranet/activation/>`_
* `SLD Calulator <https://sld-calculator.appspot.com/>`_

Documentation is `available online <https://periodictable.readthedocs.io>`_.

Source links:

* https://github.com/python-periodictable/periodictable
* https://github.com/python-periodictable/activation

|CI| |RTD| |DOI|

.. |CI| image:: https://github.com/python-periodictable/periodictable/workflows/Test/badge.svg
   :alt: Build status
   :target: https://github.com/python-periodictable/periodictable/actions

.. |DOI| image:: https://zenodo.org/badge/1146700.svg
   :alt: DOI tag
   :target: https://zenodo.org/badge/latestdoi/1146700

.. |RTD| image:: https://readthedocs.org/projects/periodictable/badge/?version=latest
   :alt: Documentation status
   :target: https://periodictable.readthedocs.io/en/latest/?badge=latest
