# This program is public domain
# Author: Paul Kienzle

"""
Extensible periodic table of elements

The periodictable package contains mass for the isotopes and density for the
elements. It calculates xray and neutron scattering information for
isotopes and elements. Composite values can be calculated from
chemical formula and density.

The table is extensible. See the user manual for details.

----

Disclaimer:

This data has been compiled from a variety of sources for the user's
convenience and does not represent a critical evaluation by the authors.
While we have made efforts to verify that the values we use match
published values, the values themselves are based on measurements
whose conditions may differ from those of your experiment.

----

"""

__docformat__ = 'restructuredtext en'
__version__ = "2.0.2"

__all__ = ['elements'] # Lazy symbols and individual elements added later

import importlib

from . import core
from . import mass
from . import density

_LAZY_MODULES: list[str] = []
_LAZY_LOAD = {
    'formula': 'formulas',
    'mix_by_weight': 'formulas',
    'mix_by_volume': 'formulas',
    'neutron_sld': 'nsf',
    'neutron_scattering': 'nsf',
    'xray_sld': 'xsf',
}
def __getattr__(name: str):
    """
    Lazy loading of modules and symbols from other modules. This is
    equivalent to using "from .formulas import formula" etc in __init__
    except that the import doesn't happen until the symbol is referenced.
    Using "from periodictable import formula" will import the symbol immediately.
    "from periodictable import *" will import all symbols, including the lazy
    """
    module_name = _LAZY_LOAD.get(name, None)
    if module_name is not None:
        # Lazy symbol: fetch name from the target module
        #print(f"from {__name__}.{module_name} import {name} [lazy]")
        module = importlib.import_module(f'{__name__}.{module_name}')
        symbol = getattr(module, name)
        globals()[name] = symbol
        return symbol
    if name in _LAZY_MODULES:
        # Lazy module: just need to import it
        #print(f"import {__name__}.{name} [lazy]")
        return importlib.import_module(f'{__name__}.{name}')
    raise AttributeError(f"module '{__name__}' has not attribute '{name}'")
def __dir__():
    return __all__
# Support 'from periodictable import *' and 'dir(periodictable)'
__all__ = [*__all__, *_LAZY_MODULES, *_LAZY_LOAD.keys()]

# Always make mass and density available
elements = core.PUBLIC_TABLE
mass.init(elements)
density.init(elements)
del mass, density

# For type hinting with vscode we need an explicit list of elements.
hydrogen = H = elements.symbol("H")
helium = He = elements.symbol("He")
lithium = Li = elements.symbol("Li")
beryllium = Be = elements.symbol("Be")
boron = B = elements.symbol("B")
carbon = C = elements.symbol("C")
nitrogen = N = elements.symbol("N")
oxygen = O = elements.symbol("O")
fluorine = F = elements.symbol("F")
neon = Ne = elements.symbol("Ne")
sodium = Na = elements.symbol("Na")
magnesium = Mg = elements.symbol("Mg")
aluminum = Al = elements.symbol("Al")
silicon = Si = elements.symbol("Si")
phosphorus = P = elements.symbol("P")
sulfur = S = elements.symbol("S")
chlorine = Cl = elements.symbol("Cl")
argon = Ar = elements.symbol("Ar")
potassium = K = elements.symbol("K")
calcium = Ca = elements.symbol("Ca")
scandium = Sc = elements.symbol("Sc")
titanium = Ti = elements.symbol("Ti")
vanadium = V = elements.symbol("V")
chromium = Cr = elements.symbol("Cr")
manganese = Mn = elements.symbol("Mn")
iron = Fe = elements.symbol("Fe")
cobalt = Co = elements.symbol("Co")
nickel = Ni = elements.symbol("Ni")
copper = Cu = elements.symbol("Cu")
zinc = Zn = elements.symbol("Zn")
gallium = Ga = elements.symbol("Ga")
germanium = Ge = elements.symbol("Ge")
arsenic = As = elements.symbol("As")
selenium = Se = elements.symbol("Se")
bromine = Br = elements.symbol("Br")
krypton = Kr = elements.symbol("Kr")
rubidium = Rb = elements.symbol("Rb")
strontium = Sr = elements.symbol("Sr")
yttrium = Y = elements.symbol("Y")
zirconium = Zr = elements.symbol("Zr")
niobium = Nb = elements.symbol("Nb")
molybdenum = Mo = elements.symbol("Mo")
technetium = Tc = elements.symbol("Tc")
ruthenium = Ru = elements.symbol("Ru")
rhodium = Rh = elements.symbol("Rh")
palladium = Pd = elements.symbol("Pd")
silver = Ag = elements.symbol("Ag")
cadmium = Cd = elements.symbol("Cd")
indium = In = elements.symbol("In")
tin = Sn = elements.symbol("Sn")
antimony = Sb = elements.symbol("Sb")
tellurium = Te = elements.symbol("Te")
iodine = I = elements.symbol("I")
xenon = Xe = elements.symbol("Xe")
cesium = Cs = elements.symbol("Cs")
barium = Ba = elements.symbol("Ba")
lanthanum = La = elements.symbol("La")
cerium = Ce = elements.symbol("Ce")
praseodymium = Pr = elements.symbol("Pr")
neodymium = Nd = elements.symbol("Nd")
promethium = Pm = elements.symbol("Pm")
samarium = Sm = elements.symbol("Sm")
europium = Eu = elements.symbol("Eu")
gadolinium = Gd = elements.symbol("Gd")
terbium = Tb = elements.symbol("Tb")
dysprosium = Dy = elements.symbol("Dy")
holmium = Ho = elements.symbol("Ho")
erbium = Er = elements.symbol("Er")
thulium = Tm = elements.symbol("Tm")
ytterbium = Yb = elements.symbol("Yb")
lutetium = Lu = elements.symbol("Lu")
hafnium = Hf = elements.symbol("Hf")
tantalum = Ta = elements.symbol("Ta")
tungsten = W = elements.symbol("W")
rhenium = Re = elements.symbol("Re")
osmium = Os = elements.symbol("Os")
iridium = Ir = elements.symbol("Ir")
platinum = Pt = elements.symbol("Pt")
gold = Au = elements.symbol("Au")
mercury = Hg = elements.symbol("Hg")
thallium = Tl = elements.symbol("Tl")
lead = Pb = elements.symbol("Pb")
bismuth = Bi = elements.symbol("Bi")
polonium = Po = elements.symbol("Po")
astatine = At = elements.symbol("At")
radon = Rn = elements.symbol("Rn")
francium = Fr = elements.symbol("Fr")
radium = Ra = elements.symbol("Ra")
actinium = Ac = elements.symbol("Ac")
thorium = Th = elements.symbol("Th")
protactinium = Pa = elements.symbol("Pa")
uranium = U = elements.symbol("U")
neptunium = Np = elements.symbol("Np")
plutonium = Pu = elements.symbol("Pu")
americium = Am = elements.symbol("Am")
curium = Cm = elements.symbol("Cm")
berkelium = Bk = elements.symbol("Bk")
californium = Cf = elements.symbol("Cf")
einsteinium = Es = elements.symbol("Es")
fermium = Fm = elements.symbol("Fm")
mendelevium = Md = elements.symbol("Md")
nobelium = No = elements.symbol("No")
lawrencium = Lr = elements.symbol("Lr")
rutherfordium = Rf = elements.symbol("Rf")
dubnium = Db = elements.symbol("Db")
seaborgium = Sg = elements.symbol("Sg")
bohrium = Bh = elements.symbol("Bh")
hassium = Hs = elements.symbol("Hs")
meitnerium = Mt = elements.symbol("Mt")
darmstadtium = Ds = elements.symbol("Ds")
roentgenium = Rg = elements.symbol("Rg")
copernicium = Cn = elements.symbol("Cn")
nihonium = Nh = elements.symbol("Nh")
flerovium = Fl = elements.symbol("Fl")
moscovium = Mc = elements.symbol("Mc")
livermorium = Lv = elements.symbol("Lv")
tennessine = Ts = elements.symbol("Ts")
oganesson = Og = elements.symbol("Og")
deuterium = D = elements.symbol("D")
tritium = T = elements.symbol("T")
neutron = n = elements.symbol("n")

# Add element name and symbol (e.g. nickel and Ni) to the public attributes.
__all__ += core.define_elements(elements, globals())

# Lazy loading of element and isotope attributes, e.g., Ni.covalent_radius
def _load_covalent_radius():
    """
    covalent radius: average atomic radius when bonded to C, N or O.
    """
    from . import covalent_radius
    covalent_radius.init(elements)
core.delayed_load(['covalent_radius',
                   'covalent_radius_units',
                   'covalent_radius_uncertainty'],
                  _load_covalent_radius)
del _load_covalent_radius

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.

    Reference:
        *Ashcroft and Mermin.*
    """
    from . import crystal_structure
    crystal_structure.init(elements)
core.delayed_load(['crystal_structure'], _load_crystal_structure)
del _load_crystal_structure

def _load_neutron():
    """
    Neutron scattering factors, *nuclear_spin* and *abundance*
    properties for elements and isotopes.

    Reference:
        *Rauch. H. and Waschkowski. W., ILL Nuetron Data Booklet.*
    """
    from . import nsf
    nsf.init(elements)
core.delayed_load(['neutron'], _load_neutron, isotope=True)
del _load_neutron

def _load_neutron_activation():
    """
    Neutron activation calculations for isotopes and formulas.

    Reference:
        *IAEA 273: Handbook on Nuclear Activation Data.*
        *NBSIR 85-3151: Compendium of Benchmark Neutron Field.*
    """
    from . import activation
    activation.init(elements)
core.delayed_load(['neutron_activation'], _load_neutron_activation,
                  element=False, isotope=True)
del _load_neutron_activation

def _load_xray():
    """
    X-ray scattering properties for the elements.

    Reference:
        *Center for X-Ray optics. Henke. L., Gullikson. E. M., and Davis. J. C.*
    """
    from . import xsf
    xsf.init(elements)
core.delayed_load(['xray'], _load_xray, ion=True)
del _load_xray

def _load_emission_lines():
    """
    X-ray emission lines for various elements, including Ag, Pd, Rh, Mo,
    Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti. *K_alpha* is the average of
    K_alpha1 and K_alpha2 lines.
    """
    from . import xsf
    xsf.init_spectral_lines(elements)
core.delayed_load(['K_alpha', 'K_beta1', 'K_alpha_units', 'K_beta1_units'],
                  _load_emission_lines)
del _load_emission_lines

def _load_magnetic_ff():
    """
    Magnetic Form Fators. These values are directly from CrysFML.

    Reference:
        *Brown. P. J.(Section 4.4.5)
        International Tables for Crystallography Volume C, Wilson. A.J.C.(ed).*
    """
    from . import magnetic_ff
    magnetic_ff.init(elements)
core.delayed_load(['magnetic_ff'], _load_magnetic_ff)
del _load_magnetic_ff


# Data needed for setup.py when bundling the package into an exe
def data_files():
    """
    Return the data files associated with all periodic table attributes.

    The format is a list of (directory, [files...]) pairs which can be
    used directly in setup(..., data_files=...) for setup.py.
    """
    import os
    import glob
    from .core import get_data_path

    def _finddata(ext, patterns):
        files = []
        path = get_data_path(ext)
        for p in patterns:
            files += glob.glob(os.path.join(path, p))
        return files

    files = [('periodictable-data/xsf',
              _finddata('xsf', ['*.nff', 'read.me'])),
             ('periodictable-data', _finddata('.', ['activation.dat', 'f0_WaasKirf.dat']))]
    return files
