# This program is in the public domain
'''
Crystal structure data.

Adds *crystal_structure* to the periodic table.  Each crystal structure
is a dictionary which contains the key 'symmetry'.  Depending on the
value of crystal_structure['symmetry'], one or more parameters
'a', 'c/a', 'b/a', 'd', and 'alpha' may be present according to
the following table:

.. table:: Crystal lattice parameters

    ============ ===========
    Symmetry     Parameters
    ============ ===========
    atom
    diatom       d
    BCC          a
    fcc          a
    hcp          c/a, a
    Tetragonal   c/a, a
    Cubic        a
    Diamond      a
    Orthorhombic c/a, a, b/a
    Rhombohedral a, alpha
    SC           a
    Monoclinic
    ============ ===========

Example:

.. doctest::

    >>> import periodictable as elements
    >>> print(elements.C.crystal_structure['symmetry'])
    Diamond
    >>> print(elements.C.crystal_structure['a'])
    3.57

This data is from Ashcroft and Mermin.
'''

from typing import Any
from .core import PeriodicTable

crystal_structures: list[dict[str, Any]|None] = [
    {'symmetry': 'diatom', 'd': 0.74}, #H
    {'symmetry': 'atom'}, #He
    {'symmetry': 'BCC', 'a': 3.49}, #Li
    {'symmetry': 'hcp', 'c/a': 1.567, 'a': 2.29}, #Be
    {'symmetry': 'Tetragonal', 'c/a': 0.576, 'a': 8.73}, #B
    {'symmetry': 'Diamond', 'a': 3.57}, #C
    {'symmetry': 'diatom', 'd': 1.10}, #N
    {'symmetry': 'diatom', 'd': 1.21}, #O
    {'symmetry': 'diatom', 'd': 1.42}, #F
    {'symmetry': 'fcc', 'a': 4.43}, #Ne
    {'symmetry': 'BCC', 'a': 4.23}, #Na
    {'symmetry': 'hcp', 'c/a': 1.624, 'a': 3.21}, #Mg
    {'symmetry': 'fcc', 'a': 4.05}, #Al
    {'symmetry': 'Diamond', 'a': 5.43}, #Si
    {'symmetry': 'Cubic', 'a': 7.17}, #P
    {'symmetry': 'Orthorhombic', 'c/a': 2.339, 'a': 10.47, 'b/a': 1.229}, #S
    {'symmetry': 'Orthorhombic', 'c/a': 1.324, 'a': 6.24, 'b/a': 0.718}, #Cl
    {'symmetry': 'fcc', 'a': 5.26}, #Ar
    {'symmetry': 'BCC', 'a': 5.23}, #K
    {'symmetry': 'fcc', 'a': 5.58}, #Ca
    {'symmetry': 'hcp', 'c/a': 1.594, 'a': 3.31}, #Sc
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 2.95}, #Ti
    {'symmetry': 'BCC', 'a': 3.02}, #V
    {'symmetry': 'BCC', 'a': 2.88}, #Cr
    {'symmetry': 'Cubic', 'a': 8.89}, #Mn
    {'symmetry': 'BCC', 'a': 2.87}, #Fe
    {'symmetry': 'hcp', 'c/a': 1.622, 'a': 2.51}, #Co
    {'symmetry': 'fcc', 'a': 3.52}, #Ni
    {'symmetry': 'fcc', 'a': 3.61}, #Cu
    {'symmetry': 'hcp', 'c/a': 1.856, 'a': 2.66}, #Zn
    {'symmetry': 'Orthorhombic', 'c/a': 1.695, 'a': 4.51, 'b/a': 1.001}, #Ga
    {'symmetry': 'Diamond', 'a': 5.66}, #Ge
    {'symmetry': 'Rhombohedral', 'a': 4.13, 'alpha': 54.10}, #As
    {'symmetry': 'hcp', 'c/a': 1.136, 'a': 4.36}, #Se
    {'symmetry': 'Orthorhombic', 'c/a': 1.307, 'a': 6.67, 'b/a': 0.672}, #Br
    {'symmetry': 'fcc', 'a': 5.72}, #Kr
    {'symmetry': 'BCC', 'a': 5.59}, #Rb
    {'symmetry': 'fcc', 'a': 6.08}, #Sr
    {'symmetry': 'hcp', 'c/a': 1.571, 'a': 3.65}, #Y
    {'symmetry': 'hcp', 'c/a': 1.593, 'a': 3.23}, #Zr
    {'symmetry': 'BCC', 'a': 3.30}, #Nb
    {'symmetry': 'BCC', 'a': 3.15}, #Mo
    {'symmetry': 'hcp', 'c/a': 1.604, 'a': 2.74}, #Tc
    {'symmetry': 'hcp', 'c/a': 1.584, 'a': 2.70}, #Ru
    {'symmetry': 'fcc', 'a': 3.80}, #Rh
    {'symmetry': 'fcc', 'a': 3.89}, #Pd
    {'symmetry': 'fcc', 'a': 4.09}, #Ag
    {'symmetry': 'hcp', 'c/a': 1.886, 'a': 2.98}, #Cd
    {'symmetry': 'Tetragonal', 'c/a': 1.076, 'a': 4.59}, #In
    {'symmetry': 'Tetragonal', 'c/a': 0.546, 'a': 5.82}, #Sn
    {'symmetry': 'Rhombohedral', 'a': 4.51, 'alpha': 57.60}, #Sb
    {'symmetry': 'hcp', 'c/a': 1.330, 'a': 4.45}, #Te
    {'symmetry': 'Orthorhombic', 'c/a': 1.347, 'a': 7.27, 'b/a': 0.659}, #I
    {'symmetry': 'fcc', 'a': 6.20}, #Xe
    {'symmetry': 'BCC', 'a': 6.05}, #Cs
    {'symmetry': 'BCC', 'a': 5.02}, #Ba
    {'symmetry': 'hcp', 'c/a': 1.619, 'a': 3.75}, #La
    {'symmetry': 'fcc', 'a': 5.16}, #Ce
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.67}, #Pr
    {'symmetry': 'hcp', 'c/a': 1.614, 'a': 3.66}, #Nd
    None, #Pm
    {'symmetry': 'Rhombohedral', 'a': 9.00, 'alpha': 23.13}, #Sm
    {'symmetry': 'BCC', 'a': 4.61}, #Eu
    {'symmetry': 'hcp', 'c/a': 1.588, 'a': 3.64}, #Gd
    {'symmetry': 'hcp', 'c/a': 1.581, 'a': 3.60}, #Tb
    {'symmetry': 'hcp', 'c/a': 1.573, 'a': 3.59}, #Dy
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.58}, #Ho
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.56}, #Er
    {'symmetry': 'hcp', 'c/a': 1.570, 'a': 3.54}, #Tm
    {'symmetry': 'fcc', 'a': 5.49}, #Yb
    {'symmetry': 'hcp', 'c/a': 1.585, 'a': 3.51}, #Lu
    {'symmetry': 'hcp', 'c/a': 1.582, 'a': 3.20}, #Hf
    {'symmetry': 'BCC', 'a': 3.31}, #Ta
    {'symmetry': 'BCC', 'a': 3.16}, #W
    {'symmetry': 'hcp', 'c/a': 1.615, 'a': 2.76}, #Re
    {'symmetry': 'hcp', 'c/a': 1.579, 'a': 2.74}, #Os
    {'symmetry': 'fcc', 'a': 3.84}, #Ir
    {'symmetry': 'fcc', 'a': 3.92}, #Pt
    {'symmetry': 'fcc', 'a': 4.08}, #Au
    {'symmetry': 'Rhombohedral', 'a': 2.99, 'alpha': 70.45}, #Hg
    {'symmetry': 'hcp', 'c/a': 1.599, 'a': 3.46}, #Tl
    {'symmetry': 'fcc', 'a': 4.95}, #Pb
    {'symmetry': 'Rhombohedral', 'a': 4.75, 'alpha': 57.14}, #Bi
    {'symmetry': 'SC', 'a': 3.35}, #Po
    None, #At
    None, #Rn
    None, #Fr
    None, #Ra
    {'symmetry': 'fcc', 'a': 5.31}, #Ac
    {'symmetry': 'fcc', 'a': 5.08}, #Th
    {'symmetry': 'Tetragonal', 'c/a': 0.825, 'a': 3.92}, #Pa
    {'symmetry': 'Orthorhombic', 'c/a': 2.056, 'a': 2.85, 'b/a': 1.736}, #U
    {'symmetry': 'Orthorhombic', 'c/a': 1.411, 'a': 4.72, 'b/a': 1.035}, #Np
    {'symmetry': 'Monoclinic'}, #Pu
    None, #Am
    None, #Cm
    None, #Bk
    None, #Cf
    None, #Es
    None, #Fm
    None, #Md
    None, #No
    None, #Lw
]

def init(table: PeriodicTable, reload: bool=False) -> None:
    """
    Add crystal_structure field to the element properties.
    """
    if 'crystal_structure' in table.properties and not reload:
        return
    table.properties.append('crystal_structure')

    for k, struct in enumerate(crystal_structures):
        table[k+1].crystal_structure = struct
