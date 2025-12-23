Known issues
============

* Incoherent scattering computed for contrast matched mixture in D2O_sld(),
  differs from the value that would be computed for a compound with the same
  isotope proportions and density computed in neutron_sld(). This may change
  in a future release.

* Incoherent scattering calculations for energy-dependent rare earth elements
  is underestimated. The calculation requires bound incoherent scattering
  length (b_i) but only the bound coherent scattering length (b_c) is
  included.

* For very high fluences, e.g., more than 1E16 n/cm2, the activation
  equations give erroneous results because of the precision limitations.
  If there is doubt simply do the calculation at a lower flux and
  proportion the result. This will not work for the cascade reactions,
  i.e., two neutron additions.

* Post-exposure decay for b reaction daughter nuclides do not include
  activity from the decaying parent nuclide.

Change history
==============

2026-01-xx R2.1.0
-----------------

Modified:

* Use NUBASE2020 for half-life values in activation. The following change by 20% or more:
  32Si, 38Clm, 41Ca, 79Se, 97Tc, 108Ag, 121Snm, 126Sn, 135Cs, 157Tb, 177Lum

2024-12-19 R2.0.2
-----------------

Modified:

* Revert bound coherent scattering length updates for C-12, O-16, Sn-119,
  Pb-207 from R2.0.0 and instead update nat-C, nat-O, nat-Sn, and nat-Pb.

2024-12-06 R2.0.1
-----------------

* Documentation tweaks

2024-12-03 R2.0.0
-----------------

Modified:

* Move to IAEA AME2020 for isotope mass
* Move to IUPAC CIAAW 2021 for atomic mass and isotopic abundance
* Move old mass/abundance tables into periodictable.mass_2001
* Li-6:Li-7 mass ratio changed from 12.2 to 19.6 (delta = 2.7%)
* Isotope abundance changed by 0.1 to 0.5 for B, Zn, Ge, Se, Mo, Er, Yb, Pt, Hg
* Atomic weight changed by 0.04% for Zn, 0.02% for S and 0.01% for Li, Ge, Se, Mo
* Fix typos in uncertainties in the neutron table (Zr-90, Te-124, Ba-138, Sm-147)
* Use correct activation for Pb-208 which was entered as barns rather than mbarns
* Use Aug 2023 update to X-ray scattering factors (Pt, Cr, Nb, Y, Er)
* Update physical constants to CODATA 2022. Values have changed by as much
  as 1e-7 (relative) in Avogadro number, Planck constant, J/eV, electron radius,
  electron mass, neutron mass and atomic mass constant.
* Updated bound coherent scattering length for H-1, H-2, He-4, C-12†,
  O-16†, O-17, O-18, Sn-119†, Sm-154, Eu-153, Ir-191, Ir-193, Pb-207†, Bi-209
* Neutron b_c changed for Zn-70 from 6.9 to 6.0 (fixes a typo in the table)
* Updated total cross section for He, Kr, Xe
* Use log-log interpolation for X-ray f''

† Some scattering length updates were reverted in R2.0.2. They should have
been updates to the element in natural abundance rather than the isotope.

Breaking changes:

* Remove Z=0 (bare neutron) from the PeriodicTable iterator.
* Replace `constants.plancks_constant` (eV s) with
  `constants.planck_constant` (J/Hz).

2024-07-08 R1.7.1
-----------------

Updated:

* Support numpy 2.x, which removed the alias np.NaN

2024-03-22 R1.7.0
-----------------

New:

* Use wt% and vol% for mixtures. The non-standard %wt and %vol are still
  supported, but may be removed in a future version if they cause ambiguity
  in the parser.
* Support mixtures containing FASTA components such as dna:CGCTAATC
* Support unicode subscripts in chemical formula input.

Modified:

* Fasta calculations for formula, density and sld have changed. DNA/RNA now
  use Buckin (1989) for unit volumes. DNA/RNA no longer include sodium ions
  from tables in Perkins (1988). Sequences now include  H+ and OH- terminators.
* Use correct halflife for Tm-171, Ho-163 and W-188 activation products.
* Fix decay time estimation routine.
* Tested on Python 3.8 and above. Support for python 2.7 dropped.
* Remove eval() from codebase.

2022-05-18 R1.6.1
-----------------

Modified:

* Calculate decay time correctly in the presence of significant long-lived
  activation.
* Calculate b mode activation correctly for Eu[151] => Gd[152]

2021-04-21 R1.6.0
-----------------

New:

* Add energy dependence for rare earths (Lynn and Seeger, 1990).

Modified:

* Use complex b_c when computing the coherent cross section, leading to
  correct values of sigma_c and sigma_i for materials with large absorption.
  With this change the tabulated values for B[10] are now shown to be
  self-consistent within a few percent.

Breaking changes:

* Neutron scattering factors are returned with one value for each wavelength
  even for energy independent elements. Previous versions returned a scalar
  if the returned value was identical for each wavelength.

2020-11-04 R1.5.3
-----------------

Breaking changes:

* Fix calculation of contrast match points for biomolecules. The old
  formula used the density of H2O for the D2O sld calculation.
* Modify biomolecule support to use H[1] rather than T for labile hydrogen.
  This will result in less error when the labile formula is used in lieu
  of the natural formula or the contrast-matched formula, and make it more
  obvious from glancing at the formula that labile hydrogen is present.
* Modify *fasta.Molecule* attributes, dropping *Hmass* and *Hsld*. *Hnatural*
  has been moved to *natural_formula*. The formula with labile hydrogen is
  stored in *labile_formula*, as well as *formula* as before.

New:

* Add *replace()* method to formula to allow isotope substitution.
* Add *nsf.D2O_match()* and *nsf.D2O_sld()* functions.

Modified:

* Neutron wavelength now defaults to 1.798 A when wavelength and energy are
  both None in *neutron_sld()* and *neutron_scattering()* rather than
  throwing an assertion error.
* *table* can be passed to neutron sld calculators as the source of isotope
  information when parsing the chemical formula.
* Switch unit test framework from nose to pytest.
* Update docs.

2019-11-19 R1.5.2
-----------------

Modified:

* Carbon density changed from 2.1 to 2.2 to match CXRO, CRC and RSC. The NIST
  X-ray attenuation tables use 2.26; the Handbook of Mineralogy has 2.09-2.23.
  The Neutron Data Booklet gave the value as 1.9-2.3, and 2.1 was chosen
  from this range.  The remaining density will continue to use values from the
  Neutron Data Booklet, which cites CRC as the primary source.
* Updated references.

2019-09-09 R1.5.1
-----------------

Modified:

* fasta uses natural abundance of H for biomolecule when computing the
  D2O contrast match rather than the biomolecule with pure H[1].
* remove half-life units from column header in activation table since
  each row gives its own units.

2017-05-11 R1.5.0
-----------------

New:

* mixture by mass and volume, e.g., 5 g NaCl // 50 mL H2O@1
* multilayer materials, e.g., 5 um Si // 3 nm Cr // 8 nm Au
* add support for bio molecules with labile hydrogens
* update list of possible oxidation states to include rare states

Modified:

* fixed computation of incoherent cross section so it is consistent with
  coherent cross section and total cross section

2014-02-04 R1.4.1
-----------------

Modified:

* default density is now the isotopic density rather than the natural density

2013-12-20 v1.4.0
-----------------

* support python 3.3

2013-10-25 R1.3.10
------------------

Modified:

* fix activation calculation to ignore fast neutrons in thermal environment
* add emission spectra for remaining elements above neon

2013-04-23 R1.3.9
-----------------

Modified:

* Update requirements to pyparsing<2.0.0 (we don't support python 3 yet)

2013-04-08 R1.3.8
-----------------

New:

* formula parser supports density spec and mix by weight/mix by volume

Modified:

* py2exe/py2app wrapping now includes missing activation.dat
* skipping bad 1.3.7 build which didn't include all changes

2013-03-05 R1.3.6
-----------------

New:

* add activation decay time to neutron activation calculator

Modified:

* Change neutron scattering calculations for incoherent cross section
  to be the linear combination of the incoherent cross sections of the
  individual atoms rather than total cross section minus the coherent
  cross section.  Penetration depth of the unscattered beam still uses
  the total cross section plus the absorption cross section.

2013-02-26 R1.3.5
-----------------

New:

* formulas now report charge and mass_fraction
* formula parser accepts ions as Yy{#+} or Yy[#]{#+} for isotopes
* support neutron activation calculations
* support xray refraction index and mirror reflectivity

Modified:

* update X-ray scattering tables for Zr
* adjust ion mass for number of electrons
* ions now display as Yy{#+} rather than Yy^{#+}
* fix formula.natural_density
* fix formula.hill so C,H come first
* fix element.interatomic_distance
* formula(value=...) -> formula(compound=...)

2010-12-05 R1.3
---------------

New:

* mix_by_weight and mix_by_volume formula constructors
* use natural density to set density for isotope specific formulas
* add neutron_scattering function which returns xs, sld and penetration depth

Modified:

* need wavelength= or energy= for xray/neutron sld
* improved docs and testing

2010-04-28 R1.2
---------------

New:

* support pickle: id(H) == id(loads(dumps(H)))
* support ions, with magnetic form factors and x-ray f0 scattering factor
* support py2exe wrappers
* allow density to be calculated from structure (bcc, fcc, hcp, cubic, diamond)
* estimate molecular volume
* support private tables with some values replaced by application

Modified:

* rename package periodictable
* rename table to periodictable.elements
* neutron sld returns real and imaginary coherent and incoherent
  instead of coherent, absorption and incoherent
* bug fix: sld for H[2] was wrong when queried before sld for H.
* remove CrysFML ionic radius definitions

2009-01-20 R1.1
---------------

Modified:

* Restructure package, separating tests into different directory
* When defining table extensions, you should now do::

      from elements.core import periodic_table, Element, Isotope

  rather than::

      from elements import periodic_table
      from elements.elements import Element, Isotope
