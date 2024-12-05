.. _formula:


********************
Chemical Composition
********************

Some properties are available for groups of elements.  Groups are specified
as a chemical formula string and either density or cell volume for the crystal
structure.  While it does not provide any information about molecular
structure, a formula does provide complete control over chemical composition.

A formula string is translated into a formula using
:func:`periodictable.formulas.formula`:

* Formula strings consist of counts and atoms, where individual atoms are
  represented by periodic table symbol.  The atoms are case sensitive,
  so "CO" is different from "Co".  Here is an example of calcium carbonate:

    >>> from periodictable import formula
    >>> print(formula("CaCO3"))
    CaCO3

* Formulas can contain multiple groups separated by space or plus or by using
  parentheses.  Whole groups can have a repeat count.  The following are
  equivalent definitions of hydrated calcium carbonate:

    >>> print(formula("CaCO3+6H2O"))
    CaCO3(H2O)6
    >>> print(formula("CaCO3 6H2O"))
    CaCO3(H2O)6
    >>> print(formula("CaCO3(H2O)6"))
    CaCO3(H2O)6

* Parentheses can nest, e.g., in polyethylene glycol:

    >>> print(formula("HO ((CH2)2O)6 H"))
    HO((CH2)2O)6H

* Isotopes are represented by index, such as O[18] = :sup:`18`\ O:

    >>> print(formula("CaCO[18]3+6H2O"))
    CaCO[18]3(H2O)6

* Ions are represented by charge, such as O{2-} = O\ :sup:`2-`:

    >>> print(formula("P{5+}O{2-}4"))
    P{5+}O{2-}4

  If charge is +/- 1 then the number is optional:

    >>> print(formula("Na{+}Cl{1-}"))
    Na{+}Cl{-}

  When specifying both charge and isotope, isotope comes first:

    >>> print(formula("Fe[56]{2+}"))
    Fe[56]{2+}

  Even though the charge is on the individual atoms, the entire formula
  has a charge:

    >>> print(formula("P{5+}O{2-}4").charge)
    -3

* Counts can be integer or decimal:

    >>> print(formula("CaCO3+(3HO1.5)2"))
    CaCO3((HO1.5)3)2

* Unicode subscripts can be used for counts, with the usual decimal point:

    >>> print(formula("CaCO₃+(3HO₁.₅)₂"))
    CaCO3((HO1.5)3)2

* Print formulas with unicode using the pretty() function:

    >>> from periodictable.formulas import pretty
    >>> print(pretty(formula("CaCO3+(3HO1.5)2")))
    CaCO₃((HO₁.₅)₃)₂

* Formula density can be specified using the special '@' tag:

    >>> print(formula("NaCl@2.16").density)
    2.16

  Density gives the isotopic density of the compound, so for example, D2O
  could be specified using:

    >>> print("%.3f"%formula("D2O@1.112").density)
    1.112

  It can also be specified using the natural density of the compound,
  assuming the isotopes substitution does not change the unit cell volume:

    >>> print("%.3f"%formula("D2O@1n").density)
    1.112

  Density applies to the entire formula, so for example a D2O-H2O
  2:1 mixture (not by mass or by volume) would be:

    >>> print("%.3f"%formula("2D2O + H2O@1n").density)
    1.074

* Mass fractions use wt%, with the final portion adding to 100%:

    >>> print(formula("10wt% Fe // 15% Co // Ni"))
    FeCo1.4214Ni7.13602

  Only the first item needs to specify that it is a mass fraction,
  and the remainder can use a bare %.

* Volume fractions use vol%, with the final portion adding to 100%:

    >>> print(formula("10vol% Fe // Ni"))
    FeNi9.68121

  Only the first item needs to specify that it is a volume fraction, and
  the remainder can use a bare %.

  Volume fraction mixing is only possible if the densities are known for
  the individual components, which will require the formula density tag
  if the component is not an element.  A density estimate is given for
  the mixture but in general it will not be correct, and should be set
  explicitly for the resulting compound.

* Specific mass can be giving with count follwed by mass units:

    >>> print(formula("5g NaCl // 50mL H2O@1"))
    NaCl(H2O)32.4395

  Density will be required for materials given by volume.  Mass will be
  stored in the *total_mass* attribute of the resulting formula.

* Multilayers can be specified by thickness:

    >>> print(formula("1 um Si // 5 nm Cr // 10 nm Au"))
    Si119.992CrAu1.41722

  Density will be required for each layer. Thickness will be stored in
  the *total_thickness* attribute of the resulting formula. Thickness can
  be converted to *total_volume* by multiplying by cross section, and to
  *total_mass* by multiplying that by *density*.

* Mixtures can nest.  The following is a 10% salt solution by weight mixed
  20:80 by volume with D2O:

    >>> print(formula("20vol% (10 wt% NaCl@2.16 // H2O@1) // D2O@1n"))
    NaCl(H2O)29.1956(D2O)122.79

* Empty formulas are supported, e.g., for air or vacuum:

    >>> print(formula())
    <BLANKLINE>
    >>> formula()
    formula('')

The grammar used for parsing formula strings is the following:

::

    formula    :: compound | mixture | nothing
    mixture    :: quantity | percentage
    quantity   :: number unit part ('//' number unit part)*
    percentage :: number 'wt%|vol%' part ('//' number '%' part)* '//' part
    part       :: compound | '(' mixture ')'
    compound   :: (composite | fasta) density?
    fasta      :: ('dna' | 'rna' | 'aa') ':' [A-Z -*]+
    composite  :: group (separator group)*
    group      :: number element+ | '(' formula ')' number
    element    :: symbol isotope? ion? number?
    symbol     :: [A-Z][a-z]*
    isotope    :: '[' integer ']'
    ion        :: '{' integer? [+-] '}'
    density    :: '@' number [ni]?
    number     :: integer | fraction
    integer    :: [1-9][0-9]*
    fraction   :: ([1-9][0-9]* | 0)? '.' [0-9]*
    separator  :: space? '+'? space?
    unit       :: mass | volume | length
    mass       :: 'kg' | 'g' | 'mg' | 'ug' | 'ng'
    volume     :: 'L' | 'mL' | 'uL' | 'nL'
    length     :: 'cm' | 'mm' | 'um' | 'nm'

Formulas can also be constructed from atoms or other formulas:

* A simple formula can be created from a bare atom:

    >>> from periodictable import Ca, C, O, H
    >>> print(formula(Ca))
    Ca

* More complex structures will require a sequences of counts and fragments.
  The fragment itself can be a structure:

    >>> print(formula( [ (1,Ca), (1,C), (3,O), (6,[(2,H),(1,O)]) ] ))
    CaCO3(H2O)6

* Structures can also be built with simple formula math:

    >>> print(formula("CaCO3") + 6*formula("H2O"))
    CaCO3(H2O)6

* Formulas can be easily cloned:

    >>> print(formula(formula("CaCO3+6H2O")))
    CaCO3(H2O)6

Density
-------

Density can be specified directly when the formula is created, or updated
within a formula.  For isotope specific formulas, the density can be given
either as the density of the formula using naturally occurring abundance
if the unit cell is approximately the same, or using the density specific
to those isotopes used.

This makes heavy water density easily specified as:

    >>> D2O = formula("D2O", natural_density=1)
    >>> print(f"{D2O} {D2O.density:.4g}")
    D2O 1.112

Density can also be estimated from the volume of the unit cell, either
by using the covalent radii of the constituent atoms and assuming some
packing factor, or by knowing the lattice parameters of the crystal
which makes up the material.  Standard packing factors for hcp, fcc,
bcc, cubic and diamond on uniform spheres can be used if the components
are of about the same size.  The formula should specify the number of
atoms in the unit cell, which is 1 for cubic, 2 for bcc and 4 for fcc.
Be sure to use the molecular mass (M.molecular_mass in g) rather
than the molar mass (M.mass in u = g/mol) in your calculations.

Because the packing fraction method relies on the covalent radius
estimate it is not very accurate:

    >>> from periodictable import elements, formula
    >>> Fe = formula("2Fe")  # bcc lattice has 2 atoms per unit cell
    >>> Fe.density = Fe.molecular_mass/Fe.volume('bcc')
    >>> print(f"{Fe.density:.3g}")
    6.55
    >>> print(f"{elements.Fe.density:.3g}")
    7.87

Using lattice parameters the results are much better:

    >>> Fe.density = Fe.molecular_mass/Fe.volume(a=2.8664)
    >>> print(f"{Fe.density:.3g}")
    7.88

Mixtures
--------

Mixtures can be created by weight or volume ratios, with the density of
the result computed from the density of the materials.  For example, the
following is a 2:1 mixture of water and heavy water:

    >>> from periodictable import formula, mix_by_volume, mix_by_weight
    >>> H2O = formula('H2O',natural_density=1)
    >>> D2O = formula('D2O',natural_density=1)
    >>> mix = mix_by_volume(H2O,2,D2O,1)
    >>> print(f"{mix} {mix.density:.4g}")
    (H2O)2D2O 1.037

Note that this is different from a 2:1 mixture by weight:

    >>> mix = mix_by_weight(H2O,2,D2O,1)
    >>> print(f"{mix} {mix.density:.4g}")
    (H2O)2.22339D2O 1.035

Except in the simplest of cases, the density of the mixture cannot be
computed from the densities of the components, and the resulting density
should be set explicitly.

Derived values
--------------

Once a formula has been created, it can be used for summary calculations.
The following is an example of hydrated quartz, which shows how to
compute molar mass and neutron/xray scattering length density:

    >>> import periodictable
    >>> SiO2 = periodictable.formula('SiO2')
    >>> hydrated = SiO2 + periodictable.formula('3H2O')
    >>> print(f"{hydrated} mass {hydrated.mass:.3f}")
    SiO2(H2O)3 mass 114.128
    >>> rho,mu,inc = periodictable.neutron_sld('SiO2+3H2O',density=1.5,wavelength=4.75)
    >>> print(f"{hydrated} neutron sld {rho:.3g}")
    SiO2(H2O)3 neutron sld 0.849
    >>> rho,mu = periodictable.xray_sld(hydrated,density=1.5,
    ... wavelength=periodictable.Cu.K_alpha)
    >>> print(f"{hydrated} X-ray sld {rho:.3g}")
    SiO2(H2O)3 X-ray sld 13.5

Biomolecules
------------

The :mod:`periodictable.fasta` module can be used to load and manage bio
molecules.  These can be used to compute molecular weights, approximate
volumes and scattering for various lipids and proteins.  In addition it
supports labile hydrogen calculations, allowing you to compute the neutron
scattering length density of the molecule in the presence of D2O as a
solvent, assuming all labile hydrogens are substituted.
