# This program is public domain
# Author: Paul Kienzle
"""
Various fundamental constants.
"""

#: Avogadro's number (mol^-1) [CODATA 2022]
avogadro_number = 6.02214076e23 #(exact) 1/mol [CODATA 2022]
#: Planck constant (eV s) [CODATA 2022]
planck_constant = 6.62607015e-34 #(exact) J/Hz [CODATA 2022]
#: Electron volt (J/eV) [CODATA 2022]
electron_volt = 1.602176634e-19 #(exact) J / eV [CODATA 2022]
#planck_constant_eV = planck_constant / electron_volt
#: speed of light c (m/s) [CODATA 2022]
speed_of_light = 299792458 #(exact) m/s [CODATA 2022]
#: electron radius r_e (m) [CODATA 2022]
electron_radius = 2.8179403205e-15 #(13) m [CODATA 2022]

# [CODATA] From NIST Reference on Constants, Units, and Uncertainty
#   http://physics.nist.gov/cuu/index.html
# [AME2020] "The Ame2020 atomic mass evaluation (II)"
#     by M.Wang, W.J.Huang, F.G.Kondev, G.Audi and S.Naimi
#     Chinese Physics C45, 030003, March 2021.
#: neutron mass (u) [CODATA 2022]
neutron_mass     = 1.00866491606 #(40) u [CODATA 2022]
neutron_mass_unc = 0.00000000040
#neutron_mass    = 1.00866491590 #(47) u [AME 2020]
#neutron_mass    = 1.00866491597 #(43) u [CODATA 2006]
#neutron_mass    = 1.00866491595 #(49) u [CODATA 2018]
#: atomic mass constant (kg / u) [CODATA 2022]
atomic_mass_constant = 1.66053906892e-27 #(52) kg / u [CODATA 2022]
#: electron mass (u) [CODATA 2022]
electron_mass = 5.485799090441e-4 #(97) u [CODATA 2022]
