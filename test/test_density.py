from periodictable import elements

def test():
    assert elements.Cm.density == 13.51
    assert elements.Cm.density_caveat == "calculated"
    assert elements.Cm.density_units == "g/cm^3"
    assert elements.D.density == elements.H.density * elements.D.mass/elements.H.mass

    # check number density and unit cell width
    assert abs(elements.Ca.interatomic_distance - 3.5016640712645786) < 1e-10
    assert abs(elements.Fe.number_density - 8.491062108378546e+22) < 1e12
    assert elements.Ca.interatomic_distance_units == "angstrom"
    assert elements.Fe.number_density_units == "1/cm^3"

    #elements.list('symbol','density','density_caveat', format="%3s %10s   %s")

    #import mass
    #mass.init()
    #elements.list('symbol','interatomic_distance')

if __name__ == "__main__": test()
