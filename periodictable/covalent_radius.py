# -*- coding: latin-1 -*-
# This program is in the public domain
# Author: Paul Kienzle
"""

This module adds the following fields to the periodic table

*   covalent_radius
*   covalent_radius_uncertainty
*   covalent_radius_units = 'angstrom'

Use :func:`init` to initialize a private table.

Data is taken from Cordero et al., 2008 [#Cordero2008]_.  Bond specific values
(single, double, or triple) are available from Pyykkö et al., 2009 [#Pyykko2009]_,
but they are generally smaller.  The CRC Handbook uses the average of Cordero
and Pyykkö. Note that the combined Cordero/Pyykkö tables are included
herein as *periodictable.covalent_radius.CorderoPyykko*, but are not yet parsed.

The abstract of Cordero reads as follows:

    A new set of covalent atomic radii has been deduced from
    crystallographic data for most of the elements with atomic
    numbers up to 96.  The proposed radii show a well behaved
    periodic dependence that allows us to interpolate a few
    radii for elements for which structural data is lacking,
    notably the noble gases. The proposed set of radii therefore
    fills most of the gaps and solves some inconsistencies in
    currently used covalent radii.  The transition metal and
    lanthanide contractions as well as the differences in covalent
    atomic radii between low spin and high spin configurations in
    transition metals are illustrated by the proposed radii set.

.. Note::

        #. Values are averages only.  The particular radius can be highly
           dependent on oxidation state and chemical compound.

        #. The paper lists values for multiple spin states on select
           elements.  We are using sp3 for carbon and low spin for manganese,
           iron and cobalt.

        #. Elements with zero or one measurements of covalent radius are
           assigned an uncertainty of 0.00.  These are He, Ne, Pm, At, Rn,
           Fr, Ac, Pa.

        #. Elements above 96 are assigned a covalent radius and uncertainty
           of None.

        #. Radii are measured from bonds to C, N or O.  The choice of which
           compound was used is element dependent.  Details are available in
           the references.

.. [#Cordero2008] Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats,
       Marc Revés, Jorge Echeverría, Eduard Cremades, Flavia Barragán and
       Santiago Alvarez. Covalent radii revisited. Dalton Trans., 2008, 2832-2838.
       `doi:10.1039/b801115j <http://dx.doi.org/10.1039/b801115j>`_

.. [#Pyykko2009] P. Pyykkö and M. Atsumi.
      Molecular Double-Bond Covalent Radii for Elements Li-E112.
      Chemistry, A European Journal, 15, 2009, 12770-12779.
      `doi:10.1002/chem.200901472 <http://dx.doi.org/10.1002/chem.200901472>`_

"""

from .core import Element

def init(table, reload=False):
    """
    Add the covalent radius property to a private table.
    Use *reload = True* to replace the covalent radius property on an
    existing table.
    """
    if 'covalent_radius' in table.properties and not reload:
        return
    table.properties.append('covalent_radius')

    # neutron covalent radius of 0.20 A ? Not sure where that value came from
    # table[0].covalent_radius = 0.20
    Element.covalent_radius_units = 'angstrom'
    Element.covalent_radius = None
    Element.covalent_radius_uncertainty = None


    for line in Cordero.split('\n'):
        fields = line.split()
        if len(fields) == 3:
            fields += ['0', '0']  # Fill in uncertainties
        if fields[0] == '-':
            continue  # Skip alternate spin states
        Z = int(fields[0])
        #el = fields[1]
        r = float(fields[2])
        dr = float(fields[3])*0.01
        #n = int(fields[4])

        table[Z].covalent_radius = r
        table[Z].covalent_radius_uncertainty = dr

# Table of radii from Cordero.  Note that in cases where there are
# multiple spin states (C, Mn, Fe, Co) only the first spin state is used.
#
#Z,  Symbol, radius(A), uncertainty (0.01A), number of measurements
Cordero = """\
1    H    0.31    5    129
2    He    0.28
3    Li    1.28    7    5789
4    Be    0.96    3    310
5    B    0.84    3    1770
6    Csp3    0.76    1    10000
-    Csp2    0.73    2    10000
-    Csp    0.69    1    171
7    N    0.71    1    2200
8    O    0.66    2    10000
9    F    0.57    3    10000
10    Ne    0.58
11    Na    1.66    9    1629
12    Mg    1.41    7    3234
13    Al    1.21    4    3206
14    Si    1.11    2    10000
15    P    1.07    3    10000
16    S    1.05    3    10000
17    Cl    1.02    4    1987
18    Ar    1.06    10    9
19    K    2.03    12    435
20    Ca    1.76    10    2647
21    Sc    1.70    7    32
22    Ti    1.60    8    231
23    V    1.53    8    389
24    Cr    1.39    5    916
25    Mnl.s.    1.39    5    321
-     Mnh.s.    1.61    8    929
26    Fel.s.    1.32    3    336
-     Feh.s.    1.52    6    1540
27    Col.s.    1.26    3    5733
-     Coh.s.    1.50    7    780
28    Ni    1.24    4    1030
29    Cu    1.32    4    1149
30    Zn    1.22    4    443
31    Ga    1.22    3    1330
32    Ge    1.20    4    1013
33    As    1.19    4    2015
34    Se    1.20    4    1717
35    Br    1.20    3    2140
36    Kr    1.16    4    5
37    Rb    2.20    9    23
38    Sr    1.95    10    1500
39    Y    1.90    7    30
40    Zr    1.75    7    93
41    Nb    1.64    6    18
42    Mo    1.54    5    97
43    Tc    1.47    7    96
44    Ru    1.46    7    1032
45    Rh    1.42    7    458
46    Pd    1.39    6    1892
47    Ag    1.45    5    1728
48    Cd    1.44    9    19
49    In    1.42    5    546
50    Sn    1.39    4    2999
51    Sb    1.39    5    609
52    Te    1.38    4    692
53    I    1.39    3    451
54    Xe    1.40    9    2
55    Cs    2.44    11    24
56    Ba    2.15    11    3076
57    La    2.07    8    190
58    Ce    2.04    9    47
59    Pr    2.03    7    58
60    Nd    2.01    6    96
61    Pm    1.99
62    Sm    1.98    8    53
63    Eu    1.98    6    167
64    Gd    1.96    6    178
65    Tb    1.94    5    55
66    Dy    1.92    7    59
67    Ho    1.92    7    48
68    Er    1.89    6    66
69    Tm    1.90    10    15
70    Yb    1.87    8    122
71    Lu    1.87    8    61
72    Hf    1.75    10    53
73    Ta    1.70    8    88
74    W    1.62    7    219
75    Re    1.51    7    476
76    Os    1.44    4    99
77    Ir    1.41    6    131
78    Pt    1.36    5    1768
79    Au    1.36    6    114
80    Hg    1.32    5    137
81    Tl    1.45    7    291
82    Pb    1.46    5    112
83    Bi    1.48    4    51
84    Po    1.40    4    4
85    At    1.50
86    Rn    1.50
87    Fr    2.60
88    Ra    2.21    2    3
89    Ac    2.15    0    1
90    Th    2.06    6    11
91    Pa    2.00    0    1
92    U    1.96    7    57
93    Np    1.90    1    22
94    Pu    1.87    1    9
95    Am    1.80    6    11
96    Cm    1.69    3    16\
"""

# Cordero/Pyykko combined table.
# rC are the Cardero radii with uncertainty.
# r# are the single/double/triple bonded Pyykko radii.
# Entries starting with # are for specific species mentioned in Cardero.
# The CRC tables are obtained by averaging Cordero rC and Pyykko r1.
#Z  Symbol   rC(A)       r1(A)    r2(A)    r3(A)
CorderoPyykko = """\
1      H     0.31(5)     0.32
2      He    0.28        0.46
3      Li    1.28(7)     1.33     1.24
4      Be    0.96(3)     1.02     0.90     0.85
5      B     0.84(3)     0.85     0.78     0.73
6      C     0.76(1)     0.75     0.67     0.60
#sp3:        0.76(1)
#sp2:        0.73(2)
#sp:         0.69(1)
7      N     0.71(1)     0.71     0.60     0.54
8      O     0.66(2)     0.63     0.57     0.53
9      F     0.57(3)     0.64     0.59     0.53
10     Ne    0.58        0.67     0.96
11     Na    1.66(9)     1.55     1.60
12     Mg    1.41(7)     1.39     1.32     1.27
13     Al    1.21(4)     1.26     1.13     1.11
14     Si    1.11(2)     1.16     1.07     1.02
15     P     1.07(3)     1.11     1.02     0.94
16     S     1.05(3)     1.03     0.94     0.95
17     Cl    1.02(4)     0.99     0.95     0.93
18     Ar    1.06(10)    0.96     1.07     0.96
19     K     2.03(12)    1.96     1.93
20     Ca    1.76(10)    1.71     1.47     1.33
21     Sc    1.70(7)     1.48     1.16     1.14
22     Ti    1.60(8)     1.36     1.17     1.08
23     V     1.53(8)     1.34     1.12     1.06
24     Cr    1.39(5)     1.22     1.11     1.03
25     Mn    1.39(5)     1.19     1.05     1.03
#low spin:   1.39(5)
#high spin:  1.61(8)
26     Fe    1.32(3)     1.16     1.09     1.02
#low spin:   1.32(3)
#high spin:  1.52(6)
27     Co    1.26(3)     1.11     1.03     0.96
#low spin:   1.26(3)
#high spin:  1.50(7)
28     Ni    1.24(4)     1.10     1.01     1.01
29     Cu    1.32(4)     1.12     1.15     1.20
30     Zn    1.22(4)     1.18     1.20
31     Ga    1.22(3)     1.24     1.17     1.21
32     Ge    1.20(4)     1.21     1.11     1.14
33     As    1.19(4)     1.21     1.14     1.06
34     Se    1.20(4)     1.16     1.07     1.07
35     Br    1.20(3)     1.14     1.09     1.10
36     Kr    1.16(4)     1.17     1.21     1.08
37     Rb    2.20(9)     2.10     2.02
38     Sr    1.95(10)    1.85     1.57     1.39
39     Y     1.90(7)     1.63     1.30     1.24
40     Zr    1.75(7)     1.54     1.27     1.21
41     Nb    1.64(6)     1.47     1.25     1.16
42     Mo    1.54(5)     1.38     1.21     1.13
43     Tc    1.47(7)     1.28     1.20     1.10
44     Ru    1.46(7)     1.25     1.14     1.03
45     Rh    1.42(7)     1.25     1.10     1.06
46     Pd    1.39(6)     1.20     1.17     1.12
47     Ag    1.45(5)     1.28     1.39     1.37
48     Cd    1.44(9)     1.36     1.44
49     In    1.42(5)     1.42     1.36     1.46
50     Sn    1.39(4)     1.40     1.30     1.32
51     Sb    1.39(5)     1.40     1.33     1.27
52     Te    1.38(4)     1.36     1.28     1.21
53     I     1.39(3)     1.33     1.29     1.25
54     Xe    1.40(9)     1.31     1.35     1.22
55     Cs    2.44(11)    2.32     2.09
56     Ba    2.15(11)    1.96     1.61     1.49
57     La    2.07(8)     1.80     1.39     1.39
58     Ce    2.04(9)     1.63     1.37     1.31
59     Pr    2.03(7)     1.76     1.38     1.28
60     Nd    2.01(6)     1.74     1.37
61     Pm    1.99        1.73     1.35
62     Sm    1.98(8)     1.72     1.34
63     Eu    1.98(6)     1.68     1.34
64     Gd    1.96(6)     1.69     1.35     1.32
65     Tb    1.94(5)     1.68     1.35
66     Dy    1.92(7)     1.67     1.33
67     Ho    1.92(7)     1.66     1.33
68     Er    1.89(6)     1.65     1.33
69     Tm    1.90(10)    1.64     1.31
70     Yb    1.87(8)     1.70     1.29
71     Lu    1.87(8)     1.62     1.31     1.31
72     Hf    1.75(10)    1.52     1.28     1.22
73     Ta    1.70(8)     1.46     1.26     1.19
74     W     1.62(7)     1.37     1.20     1.15
75     Re    1.51(7)     1.31     1.19     1.10
76     Os    1.44(4)     1.29     1.16     1.09
77     Ir    1.41(6)     1.22     1.15     1.07
78     Pt    1.36(5)     1.23     1.12     1.10
79     Au    1.36(6)     1.24     1.21     1.23
80     Hg    1.32(5)     1.33     1.42
81     Tl    1.45(7)     1.44     1.42     1.50
82     Pb    1.46(5)     1.44     1.35     1.37
83     Bi    1.48(4)     1.51     1.41     1.35
84     Po    1.40(4)     1.45     1.35     1.29
85     At    1.50        1.47     1.38     1.38
86     Rn    1.50        1.42     1.45     1.33
87     Fr    2.60        2.23     2.18
88     Ra    2.21(2)     2.01     1.73     1.59
89     Ac    2.15        1.86     1.53     1.40
90     Th    2.06(6)     1.75     1.43     1.36
91     Pa    2.00        1.69     1.38     1.29
92     U     1.96(7)     1.70     1.34     1.18
93     Np    1.90(1)     1.71     1.36     1.16
94     Pu    1.87(1)     1.72     1.35
95     Am    1.80(6)     1.66     1.35
96     Cm    1.69(3)     1.66     1.36
97     Bk    -     1.66     1.39
98     Cf    -     1.68     1.40
99     Es    -     1.65     1.40
100    Fm    -     1.67
101    Md    -     1.73     1.39
102    No    -     1.76     1.59
103    Lr    -     1.61     1.41
104    Rf    -     1.57     1.40     1.31
105    Db    -     1.49     1.36     1.26
106    Sg    -     1.43     1.28     1.21
107    Bh    -     1.41     1.28     1.19
108    Hs    -     1.34     1.25     1.18
109    Mt    -     1.29     1.25     1.13
110    Ds    -     1.28     1.16     1.12
111    Rg    -     1.21     1.16     1.18
112    Cn    -     1.22     1.37     1.30
113    Uut   -     1.36
114    Uuq   -     1.43
115    Uup   -     1.62
116    Uuh   -     1.75
117    Uus   -     1.65
118    Uuo   -     1.57\
"""
