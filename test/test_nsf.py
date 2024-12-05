from math import sqrt, pi

import numpy as np

import periodictable
from periodictable import elements, formula, nsf
from periodictable.nsf import neutron_scattering, neutron_sld
from periodictable.constants import avogadro_number as N_A, neutron_mass

def test():
    H,He,D,O = elements.H,elements.He,elements.D,elements.O
    assert H.neutron.absorption == 0.3326
    assert H.neutron.total == 82.02
    assert H.neutron.incoherent == 80.26
    assert H.neutron.coherent == 1.7568
    assert elements.Ru[101].neutron.bp == None
    assert H[1].nuclear_spin == '1/2'
    assert H[2].nuclear_spin == '1'
    assert not H[6].neutron.has_sld()

    assert He[3].neutron.b_c_i == -1.48
    assert He[3].neutron.bm_i == -5.925

    Nb = elements.Nb
    assert Nb.neutron.absorption == Nb[93].neutron.absorption

    # Check that b_c values match abundance-weighted values
    # Note: Currently they do not for match within 5% for Ar,V,Sm or Gd
    for el in elements:
        if not hasattr(el,'neutron'): continue
        b_c = 0
        complete = True
        for iso in el:
            if iso.neutron != None:
                if iso.neutron.b_c == None:
                    complete = False
                else:
                    b_c += iso.neutron.b_c*iso.neutron.abundance/100.
        if complete and b_c != 0 and abs((b_c-el.neutron.b_c)/b_c) > 0.05:
            err = abs((b_c-el.neutron.b_c)/b_c)
            ## Printing suppressed for the release version
            #print("%2s %.3f % 7.3f % 7.3f"%(el.symbol,err,b_c,el.neutron.b_c))

    # Isotopic formula.
    M = formula('Si[30]O[18]2',density=2.2)
    sld,xs,depth = neutron_scattering(M,wavelength=4.75)
    sld2 = neutron_sld(M,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    #_summarize(M)
    #_summarize(formula('O2',density=1.14))
    assert abs(sld[0] - 3.33) < 0.01
    assert abs(sld[1] - 0) < 0.01
    #assert abs(xs[2] - 0.00292) < 0.00001   # TODO fix test
    assert abs(xs[1] - 0.00569) < 0.00001
    #assert abs(depth - 4.329) < 0.001       # TODO fix test

    # Cu/Mo K-alpha = 1.89e-5 + 2.45e-7i / 1.87e-5 + 5.16e-8i

    Ni,Si = elements.Ni, elements.Si

    # Make sure molecular calculation corresponds to direct calculation
    sld = neutron_sld('Si',density=Si.density,wavelength=4.75)
    sld2 = Si.neutron.sld(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    sld,_,_ = Si.neutron.scattering(wavelength=4.75)
    sld2 = Si.neutron.sld(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    sld,xs,depth = neutron_scattering('Si',density=Si.density,wavelength=4.75)
    sld2,xs2,depth2 = Si.neutron.scattering(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-10 for v,w in zip(xs,xs2))
    assert abs(depth-depth2) < 1e-14

    # incoherent cross sections for Ni[62] used to be negative
    sld,xs,depth = neutron_scattering('Ni[62]',density=Ni[62].density,
                                      wavelength=4.75)
    assert sld[2] == 0 and xs[2] == 0
    sld,xs,depth = Ni[62].neutron.scattering(wavelength=4.75)
    assert sld[2] == 0 and xs[2] == 0
    assert Ni[62].neutron.sld()[2] == 0


    # Test call from periodictable
    sld,xs,depth = periodictable.neutron_scattering('H2O',density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('H2O',density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-10 for v,w in zip(xs,xs2))
    assert depth==depth2
    sld = periodictable.neutron_sld('H2O',density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    # Check empty formula
    sld,xs,depth = neutron_scattering('',density=0,wavelength=4.75)
    assert all(v == 0 for v in sld)
    assert all(v == 0 for v in xs)
    assert np.isinf(depth)

    # Check density == 0 works
    sld,xs,depth = neutron_scattering('Si',density=0,wavelength=4.75)
    assert all(v == 0 for v in sld)
    assert all(v == 0 for v in xs)
    assert np.isinf(depth)

    # Test natural density
    D2O_density = (2*D.mass + O.mass)/(2*H.mass + O.mass)
    sld,xs,depth = neutron_scattering('D2O',natural_density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('D2O',density=D2O_density,wavelength=4.75)
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

    # Test that sld depends on density not on the size of the unit cell
    D2O_density = (2*D.mass + O.mass)/(2*H.mass + O.mass)
    sld,xs,depth = neutron_scattering('D2O',natural_density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('2D2O',natural_density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

    # Test energy <=> velocity <=> wavelength
    # PAK: value changes with updated neutron and atomic mass constants [2024-10]
    assert abs(nsf.neutron_wavelength_from_velocity(2200) - 1.7981972755018132) < 1e-14
    assert abs(nsf.neutron_wavelength(25) - 1.8) < 0.1
    assert abs(nsf.neutron_energy(nsf.neutron_wavelength(25)) - 25) < 1e-14

    # Confirm scattering functions accept energy and wavelength
    sld,xs,depth = neutron_scattering('H2O',density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('H2O',density=1,energy=nsf.neutron_energy(4.75))
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

def test_bare_neutron():
    n = elements.n
    assert n == elements[0]
    assert n == periodictable.neutron
    n_iso = elements[0][1]
    assert n.mass == neutron_mass
    assert n_iso.mass == neutron_mass
    assert n.neutron.b_c == -37.0
    assert n.density is None
    assert n.number_density is None
    assert n.neutron.scattering()[0] is None


def test_formula():
    density = 2.52
    M = formula('B4C', density=density)
    sld,xs,depth = neutron_scattering(M,wavelength=4.75)
    # Compare to Alan Munter's numbers:
    #   SLD=7.65e-6 - 2.34e-7i /A^2
    #   inc,abs XS = 0.193, 222.4 / cm
    #   1/e = 0.004483 cm
    #   Cu/Mo K-alpha = 2.02e-5 + 1.93e-8i / 2.01e-5 + 4.64e-9i
    # Using lambda=1.798 rather than 1.8
    #   abs XS => 222.6
    #   1/e => 0.004478
    assert abs(sld[0]-7.649)<0.001
    assert abs(sld[1]-0.234)<0.001
    assert abs(xs[1]-222.6)<0.1
    #assert abs(xs[2]-0.193)<0.001   # TODO: fix test
    #assert abs(depth-0.004478)<0.000001 # TODO: fix test
    # Check that sld_inc and coh_xs are consistent
    #   cell_volume = (molar_mass/density) / N_A * 1e24
    #   number_density = num_atoms/cell_volume
    #   sigma_i = inc_xs/number_density
    #   sld_inc = 10*number_density * sqrt ( 100/(4*pi) * sigma_i )
    #   sld_re = 10*number_density * b_c.real
    #   sigma_c = 4*pi/100*((sld_re - 1j*sld_im)/(10*number_density))**2
    #   coh_xs = sigma_c * number_density
    molar_mass = 4*elements.B.mass + elements.C.mass
    cell_volume = (molar_mass/density) / N_A * 1e24
    Nb = 5 / cell_volume
    sld_inc = Nb*sqrt(100/(4*pi)*xs[2]/Nb)*10
    coh_xs = Nb*4*pi/100*(abs(sld[0] - 1j*sld[1])/(10*Nb))**2
    assert abs(sld[2] - sld_inc) < 1e-14
    assert abs(xs[0] - coh_xs) < 1e-14

def test_contrast_matching():
    from periodictable import fasta

    # Test constrast match holds for varying volume fractions (no labile)
    SiO2 = formula("SiO2@2.4")
    match, sld_real = nsf.D2O_match(SiO2)
    sld_0p0 = nsf.D2O_sld(SiO2, volume_fraction=0.0, D2O_fraction=match)
    sld_0p7 = nsf.D2O_sld(SiO2, volume_fraction=0.7, D2O_fraction=match)
    sld_1p0 = nsf.D2O_sld(SiO2, volume_fraction=1.0, D2O_fraction=match)
    assert np.isclose(sld_0p0[0], sld_real, 1e-14)
    assert np.isclose(sld_0p7[0], sld_real, 1e-14)
    assert np.isclose(sld_1p0[0], sld_real, 1e-14)
    assert not np.isclose(sld_1p0[0], sld_1p0[1], 1e-14)

    mol = fasta.LIPIDS["cholesteral"].labile_formula
    match, sld_real = nsf.D2O_match(mol)
    sld_0p7 = nsf.D2O_sld(mol, volume_fraction=0.7, D2O_fraction=match)
    assert np.isclose(sld_0p7[0], sld_real, 1e-14)

    # Test that labile hydrogens are being subsituted in contrast match
    # Note that D2O mixture is formed from pure D2O and pure water with
    # natural H:D ratios.
    mol = formula("C3H4H[1]NO@1.29n") # alanine
    sld_0p7 = nsf.D2O_sld(mol, volume_fraction=1., D2O_fraction=0.7)
    sld_0p7_direct = nsf.neutron_sld("C3H4H0.3D0.7NO@1.29n")
    #print(sld_0p7)
    #print(sld_0p7_direct)
    assert np.isclose(sld_0p7[0], sld_0p7_direct[0], 1e-14)
    assert np.isclose(sld_0p7[1], sld_0p7_direct[1], 1e-14)
    # Not testing incoherent since it will differ
    #assert np.isclose(sld_0p7[2], sld_0p7_direct[2], 1e-14)


def test_composite():
    from periodictable.nsf import neutron_composite_sld
    molecule = '3HSO4+1H2O+2CCl4'
    material = [formula(s) for s in ('HSO4','H2O','CCl4')]
    weight = np.array([3, 1, 2])
    calc = neutron_composite_sld(material, wavelength=4.75)
    sld1 = calc(weight, density=1.2)
    sld2 = neutron_sld(molecule, density=1.2, wavelength=4.75)
    #print(material, sld1)
    #print(molecule, sld2)
    assert all(np.isscalar(v) for v in sld1 + sld2)
    assert all(abs(v-w)<1e-14 for v, w in zip(sld1, sld2))

    # with wavelength array
    calc = neutron_composite_sld(material, wavelength=[4.75, 5, 6])
    sld3 = calc(weight, density=1.2)
    assert all(len(v) == 3 for v in sld3)
    assert all(v == w[0] for v, w in zip(sld1, sld3))

    # with length one wavelength array
    calc = neutron_composite_sld(material, wavelength=[4.75])
    sld4 = calc(weight, density=1.2)
    assert all(len(v) == 1 for v in sld4)
    assert all(v == w for v, w in zip(sld1, sld4))

def test_wavelength_array():
    from periodictable.nsf import neutron_scattering, neutron_sld
    material = formula('CCl4@1.5867')
    # scalar
    sld, xs, penetration = neutron_scattering(material, wavelength=4.75)
    assert all(np.isscalar(v) for v in sld + xs + (penetration,))
    # length 1
    sld, xs, penetration = neutron_scattering(material, wavelength=[4.75])
    assert all(len(v) == 1 for v in sld + xs + (penetration,))
    # length 3
    sld, xs, penetration = neutron_scattering(material, wavelength=[3, 4, 5])
    assert all(len(v) == 3 for v in sld + xs + (penetration,))

    # scalar
    sld, xs, penetration = elements.Cl.neutron.scattering(wavelength=4.75)
    assert all(np.isscalar(v) for v in sld + xs + (penetration,))
    # length 1
    sld, xs, penetration = elements.Cl.neutron.scattering(wavelength=[4.75])
    assert all(len(v) == 1 for v in sld + xs + (penetration,))
    # length 3
    sld, xs, penetration = elements.Cl.neutron.scattering(wavelength=[3, 4, 5])
    assert all(len(v) == 3 for v in sld + xs + (penetration,))


def test_energy_dependent():
    from periodictable.nsf import neutron_composite_sld, neutron_wavelength
    from periodictable.constants import avogadro_number as NA

    # Use Lu natural to test composite since xs are derived from composite
    # Use abundance from mass.py: 97.41% Lu[175] + 2.59% Lu[176]
    # Note: abundance uses mole fraction. DOI:10.1351/PAC-REP-10-06-02
    Lu = elements.Lu
    Lu_175_abundance, Lu_176_abundance = 97.401, 2.599
    Lu_equiv = f"Lu[175]{Lu_175_abundance:g}+Lu[176]{Lu_176_abundance:g}"
    # Note: skipping incoherent xs in returned value

    # Multiple wavelength energy dependent
    wavelength = [1, 2, 3, 6] # pair of wavelengths
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength, natural_density=Lu.density)
    sld2 = Lu.neutron.sld(wavelength=wavelength)
    # sld elements are arrays of length 4
    assert all(len(v) == 4 for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # Length 1 wavelength energy dependent
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength[:1], natural_density=Lu.density)
    sld2 = Lu.neutron.sld(wavelength=wavelength[:1])
    # sld elements are arrays of length 1
    #print("length 1", sld1, sld2)
    assert all(len(v) == 1 for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # Scalar wavelength energy dependent
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength[0], natural_density=Lu.density)
    sld2 = Lu.neutron.sld(wavelength=wavelength[0])
    # sld elements are scalars; note no .all() on the comparison
    #print("scalar", sld1, sld2)
    assert all(np.isscalar(v) for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # Check that composite sld calculator works with energy dependence and
    # various wavelength vectors.
    materials = formula('Lu[175]'), formula('Lu[176]')
    weights = np.array((Lu_175_abundance, Lu_176_abundance))

    # Multiple wavelength
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength, density=Lu.density)
    calc = neutron_composite_sld(materials, wavelength=wavelength)
    sld2 = calc(weights, density=Lu.density)
    assert all(len(v) == 4 for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # Length 1 wavelength
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength[:1], density=Lu.density)
    calc = neutron_composite_sld(materials, wavelength=wavelength[:1])
    sld2 = calc(weights, density=Lu.density)
    assert all(len(v) == 1 for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # scalar wavelength
    sld1 = neutron_sld(Lu_equiv, wavelength=wavelength[0], density=Lu.density)
    calc = neutron_composite_sld(materials, wavelength=wavelength[0])
    sld2 = calc(weights, density=Lu.density)
    assert all(np.isscalar(v) for v in sld1 + sld2)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

    # Check against Alex Grutter spreadsheet values computed from Lynn&Seeger
    wavelength = neutron_wavelength(80) # look at 80 meV in the table
    number_density = 30.254
    sld1 = 4.1508488, 2.2448468, 0.
    # reconstruct density from the given number density
    density = elements.Gd.mass*number_density*1e21/NA
    sld2 = neutron_sld("Gd", wavelength=wavelength, density=density)
    assert (abs((sld1[0]-sld2[0])/sld1[0]) < 1e-14).all()
    assert (abs((sld1[1]-sld2[1])/sld1[1]) < 1e-14).all()

def time_composite():
    from periodictable.nsf import neutron_composite_sld
    import time
    calc = neutron_composite_sld([formula(s) for s in ('HSO4','H2O','CCl4')],
                                 wavelength=4.75)
    q = np.array([3,1,2])
    N = 1000
    bits = [formula(s) for s in ('HSO4','H2O','CCl4')]
    tic = time.time()
    for i in range(N):
        sld = calc(q,density=1.2)
    toc = time.time()-tic
    print("composite %.1f us"%(toc/N*1e6))
    tic = time.time()
    for i in range(N):
        sld = neutron_sld(q[0]*bits[0]+q[1]*bits[1]+q[2]*bits[2],
                          density=1.2,wavelength=4.75)
    toc = time.time()-tic
    print("direct %.1f us"%(toc/N*1e6))


def test_abundance():
    # Check abundance totals to 0% or 100%
    for el in elements:
        if not hasattr(el,'neutron'): continue
        abundance=0
        for iso in el:
            if iso.neutron == None: continue
            if not hasattr(iso.neutron,'abundance'):
                print("abundance missing for %s"%iso)
            if iso.neutron.abundance == None:
                print("%s abundance=None"%iso)
            else:
                abundance += iso.neutron.abundance
        # TODO: abundance tables are not very good
        assert abs(abundance-100) < 1.1 or abundance==0,\
            "total abundance for %s is %.15g%%"%(el.symbol,abundance)


def _summarize(M):
    from periodictable.nsf import neutron_sld, neutron_xs
    sld = neutron_sld(M,wavelength=4.75)
    xs = neutron_xs(M,wavelength=4.75)
    print("%s sld %s"%(M,sld))
    print("%s xs %s 1/e %s"%(M,xs,1/sum(xs)))
    #return
    for el in list(M.atoms.keys()):
        print("%s density %s"%(el,el.density))
        print("%s sld %s"%(el,el.neutron.sld(wavelength=4.75)))
        print("%s xs"%el + " %.15g %.15g %.15g"%el.neutron.xs(wavelength=4.75))
        print("%s 1/e %s"%(el,1./sum(el.neutron.xs(wavelength=4.75))))

def molecule_table():
    # Table of scattering length densities for various molecules
    print("SLDS for some molecules")
    for molecule,density in [('SiO2',2.2),('B4C',2.52)]:
        atoms = formula(molecule).atoms
        rho,mu,inc = neutron_sld(atoms,density,wavelength=4.75)
        print("%s(%g g/cm**3)  rho=%.4g mu=%.4g inc=%.4g"
              %(molecule,density,rho,mu,inc))

def show_tables():
    molecule_table()
    periodictable.nsf.sld_table(4.75)
    periodictable.nsf.energy_dependent_table()
    periodictable.nsf.total_comparison_table()
    periodictable.nsf.coherent_comparison_table()
    periodictable.nsf.incoherent_comparison_table()

    print("""\
Specific elements with b_c values different from Neutron News 1992.
This is not a complete list.""")
    for sym in ['Sc','Te','Xe','Sm','Eu','Gd','W','Au','Hg']:
        el = getattr(elements,sym)
        print("%s %s %s %s %s"%(el.symbol,el.neutron.b_c,el.neutron.coherent,
              el.neutron.incoherent,el.neutron.absorption))



if __name__ == "__main__":
    #time_composite()
    #test_contrast_matching()
    test_composite()
    test_energy_dependent()
