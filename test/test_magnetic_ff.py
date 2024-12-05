from periodictable import elements

def test():
    # add_form_factor("Form", "MFE2", (  0.026300, 34.959702,  0.366800, 15.943500,  0.618800,  5.593500, -0.011900) )
    # add_form_factor("j2", "FE2 ",( 1.6490,16.559, 1.9064, 6.133, 0.5206, 2.137, 0.0035))
    A,a,B,b,C,c,D = elements.Fe.magnetic_ff[2].j0
    assert b == 15.9435
    ion = elements.Fe.ion[2]
    assert ion.magnetic_ff[ion.charge].j0[3] == b
    assert ion.magnetic_ff[ion.charge].j2[3] == 6.133

if __name__ == "__main__": test()
