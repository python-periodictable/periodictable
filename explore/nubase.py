from typing import Dict, Tuple

# Time units from shortest to longest
# Year-based units use YEAR_SECONDS conversion
YEAR_SECONDS = 31556926  # Exactly 365.2422 days
TIME_UNITS = {
    'ys': 1e-24,       # yoctosecond
    'zs': 1e-21,       # zeptosecond
    'as': 1e-18,       # attosecond
    'fs': 1e-15,       # femtosecond
    'ps': 1e-12,       # picosecond
    'ns': 1e-9,        # nanosecond
    'us': 1e-6,        # microsecond
    'ms': 1e-3,        # millisecond
    's': 1,            # second
    'm': 60,           # minute
    'h': 3600,         # hour
    'd': 3600*24,      # day
    'y': YEAR_SECONDS, # year
    'ky': 1e3 * YEAR_SECONDS,  # kiloyear
    'My': 1e6 * YEAR_SECONDS,  # megayear
    'Gy': 1e9 * YEAR_SECONDS,  # gigayear
    'Ty': 1e12 * YEAR_SECONDS, # terayear
    'Py': 1e15 * YEAR_SECONDS, # petayear
    'Ey': 1e18 * YEAR_SECONDS, # exayear
    'Zy': 1e15 * YEAR_SECONDS, # zettayear
    'Yy': 1e21 * YEAR_SECONDS, # yottayear
}

class NUBASEParser:
    def __init__(self, filename: str):
        self.filename = filename

    def parse_half_life(self, hl_str: str, unit_str: str, unc_str: str) -> float:
        """Convert NUBASE half-life notation to seconds."""
        hl_str = hl_str.strip()
        # Handle special cases
        if hl_str == 'stbl':
            return float('inf'), 0
        if hl_str == 'p-unst':
            return float('inf'), 0
        if not hl_str:
            return float('nan'), float('nan')

        # handle >val <val ~val
        orig = hl_str
        if hl_str[0] in '<>~':
            hl_str = hl_str[1:]
        if hl_str[-1] == '#':
            hl_str = hl_str[:-1]

        # Parse numeric value
        try:
            value = float(hl_str)
        except ValueError:
            raise ValueError(f"Invalid half-life value: {orig}")

        try:
            unc = float(unc_str)
        except ValueError:
            unc = 0
            # msg = f"Invalid half-life uncertainty: {orig} ± {unc_str}"
            # print(msg)
            # raise ValueError(msg)

        # Handle unit conversion
        unit = unit_str.strip()
        factor = TIME_UNITS[unit]
        return value * factor, unc * factor

    def parse_decay_modes(self, decay_str: str) -> Dict[str, Dict[str, Tuple[float, float]]]:
        """Parse decay modes and branching ratios from column 120-209."""
        modes = {}

        # Split into individual mode entries
        parts = decay_str.split(';')

        for part in parts:
            part = part.strip()
            if not part:
                continue

            part = part.replace('#', '')
            # Split mode from percentage
            if part[-1] == '?':
                # "mode ?" => (0 +/- 0)  [unobserved but theoretically possible]
                key = part.split(' ')[0]
                modes[key] = 0., 0.
            elif '<' in part:
                # "mode<value" => (0 +/- value/100)
                key, valstr = part.split('<')
                valstr = valstr.split(' ')[0] # "mode<value unc"
                modes[key] = 0., float(valstr)*0.01
            elif '>' in part:
                # "mode>value" => (value/100 +/- value/100)
                key, valstr = part.split('>')
                valstr = valstr.split(' ')[0] # "mode>value unc"
                val = float(valstr)*0.01
                modes[key] = val, val
            else:
                # "mode=value" => (value/100 +/- 0)
                # "mode=value unc" => (value/100 +/- unc/100) after correcting for digits after decimal
                part = part.split('[')[0]  # ignore e.g., IT=100[gs=100,m=0]
                part = part.replace('~', '=')  # treat A~val as if it were A=val
                key, valstr = part.split('=')
                val, unc = valstr.split(' ') if ' ' in valstr else (valstr, '0')
                if '+' in unc:  # value +unc-unc
                    plus, minus = unc[1:].split('-')
                    unc = plus if float(plus) > float(minus) else minus
                unc_exp = len(val.split('.')[1]) if '.' in val else 0
                modes[key] = (float(val)*0.01, float(unc)*10**(-unc_exp-2))

        # TODO: perhaps delete IS "mode": stable isotope natural abundance
        return modes

    def parse_line(self, line: str) -> Dict[str, any]:
        """Parse a single line from the NUBASE file."""
        def to_float(s):
            s = s.strip()
            if s == 'non-exist':
                return float('nan')
            return float(s.replace('#', '')) if s else float('nan')

        # Extract fixed-width columns
        A = int(line[0:3].strip())                  # Mass Number
        Z_str = line[4:8].strip()                   # Atomic Number with isomer flag
        Z = int(Z_str[:-1])                         # Remove isomer flag
        isomer_flag = int(Z_str[-1])                # Store isomer flag
        element = line[11:16].strip()               # Element symbol
        isomer = line[16:17]                        # Isomer letter
        mass_excess = to_float(line[18:31])         # Mass Excess in keV
        mass_uncertainty = to_float(line[31:42])
        exc_energy = to_float(line[42:54].strip())
        exc_uncertainty = to_float(line[54:65].strip())
        orig = line[65:67].strip()                  # Origin of Excitation Energy
        isomer_unc = line[67:68].strip()            # Isomer ordering uncertainty
        isomer_inv = line[68:69].strip()            # Reversed ordering indicator
        half_life = line[69:78].strip()             # Half-life value
        unit = line[78:80].strip()                  # Half-life unit
        hl_uncertainty = line[81:88].strip()        # Half-life uncertainty
        spin_parity = line[88:102].strip()          # Spin and Parity
        ensdf_year = line[102:104].strip()          # ENSDF update year
        discovery_year = line[114:118].strip()       # Year of Discovery
        decay_modes = line[119:209].strip()         # Decay Modes

        # Process parsed data
        isomer_state = {
            'uncertain_ordering': isomer_unc == '*',
            'reversed_ordering': isomer_inv == '&'
        }

        spin_parity_info = {
            'direct_measurement': '*' in spin_parity,
            'systematic_value': '#' in spin_parity,
            'value': spin_parity.replace('*', '').replace('#', '')
        }

        half_life_s, err = self.parse_half_life(half_life, unit, hl_uncertainty)
        decay_info = self.parse_decay_modes(decay_modes)

        return {
            'A': A,
            'Z': Z,
            'isomer_flag': isomer_flag,
            'isomer': isomer,
            'element': element,
            'mass_excess_keV': mass_excess,
            'mass_uncertainty_keV': mass_uncertainty,
            'excitation_energy_keV': exc_energy,
            'excitation_uncertainty_keV': exc_uncertainty,
            'origin': orig,
            'isomer_state': isomer_state,
            'spin_parity': spin_parity_info,
            'ensdf_year': ensdf_year,
            'discovery_year': discovery_year,
            'half_life': ' '.join((half_life, unit)),
            'half_life_uncertainty': hl_uncertainty,
            'half_life_s': half_life_s,
            'half_life_s_err': err,
            'decay_modes_str': decay_modes,
            'decay_modes': decay_info,
        }

    def parse_file(self) -> Dict[Tuple[int, int], Dict[str, any]]:
        """Parse the entire NUBASE file and return a dictionary of isotopes."""
        isotopes = {}

        with open(self.filename, 'r') as file:
            # Skip header lines (lines starting with '#')
            for line in file:
                if not line.startswith('#'):
                    try:
                        isotope_data = self.parse_line(line)
                        key = (isotope_data['Z'], isotope_data['A'], isotope_data['isomer_flag'])
                        isotopes[key] = isotope_data
                    except Exception as e:
                        print(f"Error parsing line: {str(e)}\n  {line.strip()}")
                        raise

        return isotopes

def main(show_table=False):
    import sys
    from pathlib import Path
    import re
    from collections import defaultdict, Counter
    from math import log10

    import numpy as np
    from periodictable import elements
    #from uncertainties import ufloat as U

    decay_modes = Counter()
    checked = {}
    nubase_path = Path(__file__).parent / "nubase_4.mas20.txt"
    nubase_url = "https://www.anl.gov/phy/reference/nubase-2020-nubase4mas20"
    if not nubase_path.exists():
        print(f"Missing {nubase_path}. Fetch it from:\n   {nubase_url}")
        sys.exit(1)

    isotopes = NUBASEParser(nubase_path).parse_file()
    #for k, v in isotopes.items():
    #    print(f"{k}: {v['element']:>5s}{v['isomer']}: {v['half_life']:>10s} ({v['half_life_s']:10.3e} s) {v['decay_modes']}")

    if show_table:
        print("# Halflife from NUBASE2020 in seconds")
        print("# F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, and G. Audi, Chin. Phys. C45, 030001 (2021)")
        print("#   isomer: (halflife, uncertainty, formatted)")
        print("nubase2020_halflife = {")

    for el in elements:
        for iso_num in el.isotopes:
            iso = el[iso_num]
            act_list = getattr(iso, 'neutron_activation', ())
            #print(iso.__dict__)
            for act in act_list:
                if act.daughter in checked:
                    continue
                checked[act.daughter] = True
                match = re.match(r'^([A-Z][a-z]*)-(\d+)(.*)', act.daughter)
                d_symbol, d_isotope, d_isomer = match.group(1), match.group(2), match.group(3)
                code = 2 if d_isomer.startswith('m2') else 1 if d_isomer.startswith('m') else 0
                code_override = {
                    "Pb-204m": 2,
                    "Ir-194m": 2,
                    "Ta-182m+": 2,
                    "Lu-177m*": 3,  # 177Lup in NUBASE
                    "Eu-152m2+": 5,  # 152Eur in NUBASE
                    "Sb-122m+": 3,  # 122Sbp in NUBASE
                    "Ag-110m*": 2,
                    "Pd-109m+": 2,
                    "Pd-107m+": 2,
                    "Kr-83m": 2,
                    "Ge-73m": 2,
                    "Sc-46m+": 2,
                }
                code = code_override.get(act.daughter, code)
                d_el = getattr(elements, d_symbol)
                v = isotopes.get((d_el.number, int(d_isotope), code), None)
                if v is None:
                    print(f"{act.daughter} {act.Thalf_str} missing from nubase")
                    continue

                # While debugging the code_override above I suppressed the known differences
                known_diff = (
                    # 20% or more
                    'Lu-177m*', 'Tb-157', 'Cs-135', 'Sn-126', 'Sn-121m',
                    'Ag-108m*', 'Tc-97', 'Se-79', 'Ca-41', 'Cl-38m+', 'Si-32',
                    # 10%-20%
                    'Be-10', 'Rb-88', 'Mo-93', 'Cd-113', 'Sn-119m', 'Te-121',
                    'Lu-176', 'Os-191m+', 'Pb-205', 'Bi-210ms',
                    )
                #if act.daughter in known_diff: continue

                nub_s, nub_ds, act_s = v['half_life_s'], v['half_life_s_err'], act.Thalf_hrs*3600.0

                # Te-123 is listed as stable in NUBASE2020
                if not np.isfinite(nub_s):
                    # print(f"{act.isotope} => {act.daughter}: T1/2 = {act.Thalf_str} != stbl")
                    continue

                # # Ignore long and short half lives (in seconds)
                # if not (600 <= nub_s <= 10*365*24*3600):
                #     continue
                decay_modes.update(v['decay_modes'].keys())

                if show_table: # print new table of half-lives
                    digits = int(log10(nub_s)) - int(log10(nub_ds)) + 1
                    print(f"    \"{act.daughter}\": ({nub_s:.{digits}e}, {nub_ds:.1e}, \"{v['half_life']}\"),") # {U(nub_s, nub_ds):.2uS}
                    continue

                # Note: IS is isotopic abundance for stable and long-lived isotopes
                # Note: 2B+ from Gd152 has 0 probability
                #if any(mode not in {'B-', 'IT', 'B+', 'EC', 'A', 'IS', '2B+'} for mode in v['decay_modes'].keys()):
                #    print(f"** rare {d_symbol}{d_isotope}{v['isomer']} => {act.reaction}", v['decay_modes'])

                diff = abs(nub_s - act_s)/nub_s # relative diff between nubase and activation calculator
                target = nub_ds/nub_s # relative uncertainty in nub_s
                target = 0.005 # 0.5%
                target = 0.05 # 5%
                target = 0.2 # 20%
                if diff > target:
                    #print(act.daughter, d_symbol, d_isotope, d_isomer, code)
                    #print(f"{act.isotope}=>{act.daughter}: T1/2 = {act.Thalf_str} ({act.Thalf_hrs*3600.:.3e} s)")
                    #print(f"   {v['element']}{v['isomer']}[{v['isomer_flag']}]: {v['half_life']:>10s} ({v['half_life_s']:10.3e} s) {v['decay_modes']}")
                    print(f"{act.isotope:6} => {act.daughter:<9}: {int(diff*100):3d}% ΔT½ ({act.Thalf_str} != {v['half_life']} ± {nub_ds/nub_s*100:.1f}%)")

    if show_table:
        print("}")
        return

    print(f"There are {len(isotopes)} halflives in NUBASE of which {len(checked)} are used in the activation calculator.")
    print(f"Decay modes used in activation isotopes:", decay_modes)

    joined = {key[:2]:value for key, value in isotopes.items()}
    masses = {(el.number, iso): el[iso] for el in elements for iso in el.isotopes}
    print(f"There are {len(joined)} isotopes in NUBASE and {len(masses)} isotopes in periodictable")
    if all(m in joined for m in masses.keys()):
        print(f"All isotopes from periodictable are in nubase")
    if 0 and any(m not in masses for m in joined.keys()):
        result = defaultdict(list)
        for key, value in joined.keys():
            result[key].append(value)
        for el in elements:
            print(f"{el}:", " ".join(f"{iso}" if iso in el.isotopes else f"*{iso}" for iso in result[el.number]))

    # Note: decay energy is the difference in excess energy


if __name__ == "__main__":
    # Tell activation.py to use activation.dat for halflife values
    import periodictable.activation
    periodictable.activation.use_nubase2020_halflife = False
    main(show_table=False) # Show differences from activation.dat values
    #main(show_table=True) # Show halflife table to include in activation.py