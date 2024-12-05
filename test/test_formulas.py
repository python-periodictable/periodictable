from __future__ import division
from copy import deepcopy
from pickle import loads, dumps

from periodictable import Ca, C, O, H, Fe, Ni, Si, D, Na, Cl, Co, Ti, S
from periodictable import formula, mix_by_weight, mix_by_volume
from periodictable.formulas import count_elements, pretty

def check_parse_fails(s):
    try:
        formula(s)
    except Exception as exc:
        return True
    raise Exception(f'formula("{s}") should fail to parse')

def test():
    ikaite = formula()
    # Note: this should be a tuple of tuples
    ikaite.structure = ((1, Ca), (1, C), (3, O), (6, ((2, H), (1, O))))

    # Test print
    assert str(ikaite) == "CaCO3(H2O)6"

    # Test constructors
    assert ikaite == formula([(1, Ca), (1, C), (3, O), (6, [(2, H), (1, O)])])
    assert ikaite == formula(ikaite)
    assert ikaite is not formula(ikaite)
    assert ikaite.structure is formula(ikaite).structure

    # Test parsers
    assert formula("Ca") == formula([(1, Ca)])
    assert formula("Ca") == formula(Ca)
    assert formula("CaCO3") == formula([(1, Ca), (1, C), (3, O)])
    assert ikaite == formula("CaCO3+6H2O")
    assert ikaite == formula("(CaCO3+6H2O)1")
    assert ikaite == formula("CaCO3 6H2O")
    assert ikaite == formula("CaCO3(H2O)6")
    assert ikaite == formula("(CaCO3(H2O)6)1")
    assert ikaite.hill == formula("CCaO3(H2O)6").hill
    assert str(ikaite.hill) == "CH12CaO9"
    assert formula([(0.75, Fe), (0.25, Ni)]) == formula("Fe0.75Ni0.25")

    # Unicode, latex and html subscripts
    assert formula([(0.75, Fe), (0.25, Ni)]) == formula("Fe₀.₇₅Ni₀.₂₅")
    assert ikaite == formula("CaCO₃(H₂O)₆")
    assert ikaite == formula("CaCO₃6H₂O") # with subscripts we know it isn't O36
    assert pretty(ikaite, 'unicode') == "CaCO₃(H₂O)₆"
    assert pretty(ikaite, 'html') == "CaCO<sub>3</sub>(H<sub>2</sub>O)<sub>6</sub>"
    assert pretty(ikaite, 'latex') == "CaCO$_{3}$(H$_{2}$O)$_{6}$"
    # Only allow subscripts in the post position
    assert check_parse_fails("₃H₂O")
    assert check_parse_fails("H₂O@₁")
    assert check_parse_fails("₁wt% NaCl@2.3 // H₂O@1n")

    # Test composition
    #print formula("CaCO3") + 6*formula("H2O")
    assert ikaite == formula("CaCO3") + 6*formula("H2O")
    f = formula('')
    assert not (3*f).structure
    f = formula('H2O')
    assert id((1*f).structure) == id(f.structure)

    # Check atom count
    assert formula("Fe2O4+3H2O").atoms == {Fe: 2, O: 7, H: 6}

    # Check element count. The formula includes element, charged element,
    # isotope and charged isotope. The "3" in front forces recursion into a
    # formula tree.
    f = formula("3HDS{6+}O{2-}3O[16]{2-}")
    assert count_elements(f) == {S: 3, O: 12, H: 6}
    assert str(formula(count_elements(f)).hill) == "H6O12S3"
    assert count_elements(f, by_isotope=True) == {S: 3, O: 9, O[16]:3, H: 3, D: 3}

    # Check charge
    assert formula("P{5+}O{2-}4").charge == -3
    try:
        formula("P{18-}")
        raise Exception("No exception raised for invalid charge")
    except ValueError:
        pass
    assert formula("Na{+}Cl{-}").charge == 0
    Na_frac = Na.ion[1].mass/(Na.ion[1].mass+Cl.ion[-1].mass)
    assert abs(formula("Na{+}Cl{-}").mass_fraction[Na.ion[1]] - Na_frac) < 1e-14

    # Check the mass calculator
    assert formula('H2O').mass == 2*H.mass+O.mass
    assert formula("Fe2O4+3H2O").mass == 2*Fe.mass+7*O.mass+6*H.mass
    assert (formula("Fe2O[18]4+3H2O").mass
            == 2*Fe.mass+4*O[18].mass+3*O.mass+6*H.mass)

    # Check natural density support
    assert (formula('D2O', natural_density=1).density
            == (2*D.mass + O.mass)/(2*H.mass + O.mass))
    D2O = formula('D2O', natural_density=1)
    D2Os = formula('D2O')
    D2Os.natural_density = 1
    assert abs(D2O.density - D2Os.density) < 1e-14
    assert abs(D2O.natural_density - 1) < 1e-14
    assert abs(D2Os.natural_density - 1) < 1e-14

    # Test isotopes; make sure this is last since it changes ikaite!
    assert ikaite != formula("CaCO[18]3+6H2O")
    assert formula("O[18]").mass == O[18].mass

    # Check x-ray and neutron sld
    rho, mu, inc = formula('Si', Si.density).neutron_sld(wavelength=4.5)
    rhoSi, muSi, incSi = Si.neutron.sld(wavelength=4.5)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14
    assert abs(inc - incSi) < 1e-14

    rho, mu = formula('Si', Si.density).xray_sld(wavelength=1.54)
    rhoSi, muSi = Si.xray.sld(wavelength=1.54)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14

    # Check that names work
    permalloy = formula('Ni8Fe2', 8.692, name='permalloy')
    assert str(permalloy) == 'permalloy'

    # Check that get/restore state works
    assert deepcopy(permalloy).__dict__ == permalloy.__dict__

    # Check that copy constructor works
    #print permalloy.__dict__
    #print formula(permalloy).__dict__
    assert formula(permalloy).__dict__ == permalloy.__dict__
    assert formula('Si', name='Silicon').__dict__ != formula('Si').__dict__

    H2O = formula('H2O', natural_density=1)
    D2O = formula('D2O', natural_density=1)
    fm = mix_by_weight(H2O, 3, D2O, 2)
    fv = mix_by_volume(H2O, 3, D2O, 2)
    # quantity of H+D should stay in 2:1 ratio with O
    assert abs(fv.atoms[H]+fv.atoms[D] - 2*fv.atoms[O]) < 1e-14
    assert abs(fm.atoms[H]+fm.atoms[D] - 2*fm.atoms[O]) < 1e-14
    # H:D ratio should match H2O:D2O ratio when mixing by volume, but should
    # be skewed toward the lighter H when mixing by mass.
    assert abs(fv.atoms[H]/fv.atoms[D] - 1.5) < 1e-14
    assert abs(fm.atoms[H]/fm.atoms[D] - 1.5*D2O.density/H2O.density) < 1e-14
    # Mass densities should average according to H2O:D2O ratio when
    # mixing by volume but be skewed toward toward the more plentiful
    # H2O when mixing by mass
    H2O_fraction = 0.6
    assert abs(fv.density - (H2O.density*H2O_fraction + D2O.density*(1-H2O_fraction))) < 1e-14
    H2O_fraction = (3/H2O.density) / (3/H2O.density + 2/D2O.density)
    assert abs(fm.density - (H2O.density*H2O_fraction + D2O.density*(1-H2O_fraction))) < 1e-14

    # Make sure we are independent of unit cell size
    H2O = formula('3.2H2O', natural_density=1)
    D2O = formula('4.1D2O', natural_density=1)
    fm = mix_by_weight(H2O, 3, D2O, 2)
    fv = mix_by_volume(H2O, 3, D2O, 2)
    # quantity of H+D should stay in 2:1 ratio with O
    assert abs(fv.atoms[H]+fv.atoms[D] - 2*fv.atoms[O]) < 1e-14
    assert abs(fm.atoms[H]+fm.atoms[D] - 2*fm.atoms[O]) < 1e-14
    # H:D ratio should match H2O:D2O ratio when mixing by volume, but should
    # be skewed toward the lighter H when mixing by mass.
    assert abs(fv.atoms[H]/fv.atoms[D] - 1.5) < 1e-14
    assert abs(fm.atoms[H]/fm.atoms[D] - 1.5*D2O.density/H2O.density) < 1e-14
    # Mass densities should average according to H2O:D2O ratio when
    # mixing by volume but be skewed toward toward the more plentiful
    # H2O when mixing by mass
    H2O_fraction = 0.6
    assert abs(fv.density - (H2O.density*H2O_fraction + D2O.density*(1-H2O_fraction))) < 1e-14
    H2O_fraction = (3/H2O.density) / (3/H2O.density + 2/D2O.density)
    assert abs(fm.density - (H2O.density*H2O_fraction + D2O.density*(1-H2O_fraction))) < 1e-14

    # Pickle test
    assert loads(dumps(fm)) == fm
    ion = Fe[56].ion[2]
    assert id(loads(dumps(ion))) == id(ion)

    # zero quantities tests in mixtures
    f = mix_by_weight(H2O, 0, D2O, 2)
    assert f == D2O
    f = mix_by_weight(H2O, 2, D2O, 0)
    assert f == H2O
    f = mix_by_weight(H2O, 0, D2O, 0)
    assert f == formula()
    f = mix_by_volume(H2O, 0, D2O, 2)
    assert f == D2O
    f = mix_by_volume(H2O, 2, D2O, 0)
    assert f == H2O
    f = mix_by_volume(H2O, 0, D2O, 0)
    assert f == formula()

    # mix by weight with unknown component density
    # can't do mix by volume without component densities
    glass = mix_by_weight('SiO2', 75, 'Na2O', 15, 'CaO', 10, density=2.52)

    # layers and mixtures
    check_formula(formula('1mm Fe // 1mm Ni'), formula('50%vol Fe // Ni'))
    # The relative quantities change whenenver the mass is updated.
    #print(formula('2mL Co // 2mL Ti').structure)
    #print(formula('2g Co // 2g Ti').structure)
    #print(formula('5g NaCl // 50mL H2O@1').structure)
    check_formula(formula('50vol% Co // Ti'), formula('2mL Co // 2mL Ti'))
    check_formula(formula('50wt% Co // Ti'), formula('2g Co // 2g Ti'))
    check_formula(formula('2mL Co // 2mL Ti'), formula(((1.5922467977437773, Co), (1, Ti))))
    check_formula(formula('2g Co // 2g Ti'), formula(((1, Co), (1.2311862870035726, Ti))))

    check_formula(formula('5g NaCl // 50mL H2O@1'), formula('5g NaCl // 50g H2O'))
    check_formula(
        formula('5g NaCl // 50mL H2O@1'),
        formula(((1, Na), (1, Cl), (32.43950556758257, ((2, H), (1, O))))), tol=1e-5)
    assert abs(formula('1mm Fe // 1mm Ni').thickness - 0.002) < 0.002*1e014
    assert abs(formula('2g Co // 2g Ti').total_mass - 4) < 4*1e-14
    check_mass(formula('2mL Co // 2mL Ti'), mass=2*(Co.density+Ti.density))
    check_mass(
        formula("50 g (49 mL H2O@1 // 1 g NaCl) // 20 mL D2O@1n"),
        mass=50 + 20*D2O.density)
    check_mass(
        formula("50 mL (45 mL H2O@1 // 5 g NaCl)@1.0707 // 20 mL D2O@1n"),
        mass=50*1.0707 + 20*D2O.density)

    # fasta
    check_formula(formula('aa:A'), formula('C3H4H[1]3NO2'))
    check_formula(formula('aa:RELEEL'), formula('C33H42H[1]13N9O13'))
    check_formula(formula('aa:RELEEL'), formula('aa:RE-LEE L *UNUSED'))
    check_formula(
        formula('30%vol CCl4@1.2 //10% aa:RE-LE EL @1.8 // H2O@1'),
        formula('30%vol CCl4@1.2 //10% C33H42H[1]13N9O13 @1.8 // H2O@1'))

def check_mass(f1, mass, tol=1e-14):
    """Check that the total mass of f1 is as expected."""
    assert abs(f1.total_mass - mass) < mass*tol

def check_formula(f1, f2, tol=1e-14):
    """Check that the number of atoms in f1 and f2 are about equal."""
    f2_atoms = f2.atoms
    for atom, count in f1.atoms.items():
        if atom not in f2_atoms or abs(f2_atoms[atom] - count) > tol*count:
            raise RuntimeError("Formulas differ: %s and %s"%(f1, f2))

if __name__ == "__main__":
    test()
