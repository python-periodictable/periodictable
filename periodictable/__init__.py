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
__version__ = "2.0.1"

__all__ = ['elements'] # Lazy symbols and individual elements added later

import importlib

from . import core
from . import mass
from . import density

_LAZY_MODULES = []
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

def _load_crystal_structure():
    """
    Add crystal_structure property to the elements.

    Reference:
        *Ashcroft and Mermin.*
    """
    from . import crystal_structure
    crystal_structure.init(elements)
core.delayed_load(['crystal_structure'], _load_crystal_structure)

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

def _load_xray():
    """
    X-ray scattering properties for the elements.

    Reference:
        *Center for X-Ray optics. Henke. L., Gullikson. E. M., and Davis. J. C.*
    """
    from . import xsf
    xsf.init(elements)
core.delayed_load(['xray'], _load_xray, ion=True)

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


# Data needed for setup.py when bundling the package into an exe
def data_files():
    """
    Return the data files associated with all periodic table attributes.

    The format is a list of (directory, [files...]) pairs which can be
    used directly in setup(..., data_files=...) for setup.py.
    """
    import os
    import glob
    def _finddata(ext, patterns):
        files = []
        path = core.get_data_path(ext)
        for p in patterns:
            files += glob.glob(os.path.join(path, p))
        return files

    files = [('periodictable-data/xsf',
              _finddata('xsf', ['*.nff', 'read.me'])),
             ('periodictable-data', _finddata('.', ['activation.dat', 'f0_WaasKirf.dat']))]
    return files
