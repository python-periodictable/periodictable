from periodictable import H, O, Fe, helium, elements, data_files

def test():
    ## Uncomment to print the package path on the CI infrastructure
    #import periodictable
    #print(f"Path to imported periodictable in test dir is {periodictable.__file__}")
    #print(fail_test)
    # Check that we can access element properties
    assert H.name == "hydrogen"
    assert H.symbol == "H"
    assert H.number == 1
    assert helium.symbol == 'He'

    # Check that isotopes work and produce the correct strings and symbols
    O.add_isotope(18)
    assert H[2].symbol == 'D'
    assert H[3].symbol == 'T'
    assert O[18].symbol == 'O'
    assert str(H[2]) == 'D'
    assert str(H[3]) == 'T'
    assert str(O[18]) == '18-O'
    try:
        Fe[12]
    except KeyError as msg:
        assert msg.args[0] == '12 is not an isotope of Fe'

    # Check that "for el in elements" works and for iso in el works
    els = tuple(el for el in elements)
    assert els[0].number == 1
    assert els[1].number == 2
    isotopes = tuple(iso for iso in O)
    assert isotopes[0].isotope == 12  # 12 is the first oxygen isotope listed

    # Check that table lookup works and fails appropriately
    Fe.add_isotope(56)
    assert elements.symbol('Fe') == Fe
    assert elements.name('iron') == Fe
    assert elements.isotope('Fe') == Fe
    assert elements.isotope('56-Fe') == Fe[56]
    assert elements.isotope('D') == H[2]
    try:
        elements.symbol('Qu')
    except ValueError as msg:
        assert str(msg) == "unknown element Qu"
    try:
        elements.name('Qu')
    except ValueError as msg:
        assert str(msg) == "unknown element Qu"
    try:
        elements.isotope('Qu')
    except ValueError as msg:
        assert str(msg) == "unknown element Qu"
    try:
        elements.isotope('4-D')
    except ValueError as msg:
        assert str(msg) == "unknown element 4-D"

    # Check that ions work
    assert Fe.ion[2].charge == 2
    assert Fe.ions == (-4, -2, -1, 1, 2, 3, 4, 5, 6, 7)
    assert str(Fe.ion[2]) == "Fe{2+}"
    assert str(O.ion[-2]) == "O{2-}"
    try:
        Fe.ion[-3]
        raise Exception("accepts invalid ions")
    except ValueError as msg:
        assert str(msg) == "-3 is not a valid charge for Fe"

    assert data_files()[0][0] == "periodictable-data/xsf"

if __name__ == "__main__":
    test()
