# This program is public domain
# Author: Paul Kienzle
"""
Biomolecule support.

:class:`Molecule` lets you define biomolecules with labile hydrogen atoms
specified using H[1] in the chemical formula.  The biomolecule object creates
forms with natural isotope ratio, all hydrogen and all deuterium. Density can be
provided as natural density or cell volume.  A %D2O contrast match value is
computed for matching the molecule SLD in the presence of labile hydrogens.
:meth:`Molecule.D2Osld` computes the neutron SLD for the solvated molecule in a
%D2O solvent.

:class:`Sequence` lets you read amino acid and DNA/RNA sequences from FASTA
files.

Tables for common molecules are provided[1]:

    *AMINO_ACID_CODES* : amino acids indexed by FASTA code

    *RNA_CODES*, DNA_CODES* : nucleic bases indexed by FASTA code

    *RNA_BASES*, DNA_BASES* : individual nucleic acid bases

    *NUCLEIC_ACID_COMPONENTS*, *LIPIDS*, *CARBOHYDRATE_RESIDUES*

Neutron SLD for water at 20C is also provided as *H2O_SLD* and *D2O_SLD*.

For unmodified protein an H and an OH are added for terminations.

Assumes that proteins were created in an environment with the usual H/D isotope
ratio on the nonlabile hydrogen.

The value of residue volumes differs from that used by the bio scattering
calculators from ISIS and ORSO, which will lead to different values for SLD.
There are small differences for the number of hydrogen in His and Cys residues,
where one table considers them present but labile and the other considers them
absent.

DNA and RNA residues from the source[1] included sodium in the chemical formula,
but these have been removed and will not appear in the sequence. Volumes for DNA
and RNA residues come from Buckin (1989) as reported in Durchlag (1997), with
correction for phosphorylation and dehydration. The correction value of 30.39
comes from comparison of the volume given in Harroun (2006) to the volumes of
the RNA ACGU and DNA T nucleosides given in Buckin (1989) after correcting for
units. Harroun doesn't give volumes for DNA AGC nucleosides despite them being
different (especially guanosine). This code uses the values from Buckin for
these as well, rather than the RNA nucleoside values given in Harroun. Note that
the computed density for equal parts AGCT is 1.67, compared to the measured
average of 1.70 given in Arrighi (1970).

[1] Perkins, S.J. (1988). Chapter 6 X-Ray and Neutron Solution Scattering, in:
New Comprehensive Biochemistry. Elsevier, pp. 143-265.
https://doi.org/10.1016/S0167-7306(08)60575-X

[2] Buckin, V. A., B. I. Kankiya, and R. L. Kazaryan (1989). Hydration of
nucleosides in dilute aqueous solutions: ultrasonic velocity and density
measurements. Biophysical chemistry 34.3 211-223.
https://doi.org/10.1016/0301-4622(89)80060-2

[3] Durchschlag, H. and Zipper, P. (1997). Calculation of partial specific
volumes and other volumetric properties of small molecules and polymers. Journal
of Applied Chemistry 30 803-807. https://doi.org/10.1107/S0021889897003348

[4] Harroun, T.A., Wignall, G.D., Katsaras, J. (2006). Neutron Scattering for
Biology. In: Neutron Scattering in Biology. Springer, Berlin, Heidelberg.
https://doi.org/10.1007/3-540-29111-3_1

[5] Arrighi, F.E., Mandel, M., Bergendahl, J. et al. (1970). Buoyant densities
of DNA of mammals. Biochem Genet 4, 367–376. https://doi.org/10.1007/BF00485753
"""
import warnings
from pathlib import Path
# Warning: name clash with Sequence
from collections.abc import Iterator
from typing import IO

from .formulas import formula as parse_formula, Formula, FormulaInput
from .nsf import neutron_sld
from .xsf import xray_sld
from .core import default_table, Atom
from .constants import avogadro_number

# CRUFT 1.5.2: retaining fasta.isotope_substitution for compatibility
def isotope_substitution(formula: Formula, source: Atom, target: Atom, portion: float=1):
    """
    Substitute one atom/isotope in a formula with another in some proportion.

    *formula* is the formula being updated.

    *source* is the isotope/element to be substituted.

    *target* is the replacement isotope/element.

    *portion* is the proportion of source which is substituted for target.

    .. deprecated:: 1.5.3
        Use formula.replace(source, target, portion) instead.
    """
    return formula.replace(source, target, portion=portion)

# TODO: allow Molecule to be used as compound in formulas.formula()
class Molecule:
    """
    Specify a biomolecule by name, chemical formula, cell volume and charge.

    Labile hydrogen positions should be coded using H[1] rather than H.
    H[1] will be substituded with H for solutions with natural water
    or D for solutions with heavy water. Any deuterated non-labile hydrogen
    can be marked with D, and they will stay as D regardless of the solvent.

    *name* is the molecule name.

    *formula* is the chemical formula as string or atom dictionary, with
    H[1] for labile hydrogen.

    *cell_volume* is the volume of the molecule. If None, cell volume will
    be inferred from the natural density of the molecule. Cell volume is
    assumed to be independent of isotope.

    *density* is the natural density of the molecule. If None, density will
    be inferred from cell volume.

    *charge* is the overall charge on the molecule. Note that charge can be
    fractional if the molecule is the average of a statistical ensemble.

    **Attributes**

    *labile_formula* is the original formula, with H[1] for the labile H.
    You can retrieve the deuterated from using::

        molecule.labile_formula.replace(elements.H[1], elements.D)

    *natural_formula* has H substituted for H[1] in *labile_formula*.

    *D2Omatch* is percentage of D2O by volume in H2O required to match the
    SLD of the molecule, including substitution of labile hydrogen in
    proportion to the D/H ratio in the solvent. Values will be outside
    the range [0, 100] if the contrast match is impossible.

    *sld*/*Dsld* are the the SLDs of the molecule with H[1] replaced by
    naturally occurring H/D ratios and pure D respectively.

    *mass*/*Dmass* are the masses for natural H/D and pure D respectively.

    *charge* is the charge on the molecule

    *cell_volume* is the estimated cell volume for the molecule

    *density* is the estimated molecule density

    Change 1.5.3: drop *Hmass* and *Hsld*. Move *formula* to *labile_formula*.
    Move *Hnatural* to *formula*.
    """
    name: str
    cell_volume: float
    sld: float
    Dsld: float
    mass: float
    Dmass: float
    D2Omatch: float
    charge: float # fractional to allow ensemble average molecules
    natural_formula: Formula
    labile_formula: Formula
    formula: Formula
    code: str  # fasta code for dna/rna

    def __init__(
            self,
            name: str,
            formula: FormulaInput,
            cell_volume: float|None=None,
            density: float|None=None,
            charge: float=0,
            ):
        # TODO: fasta does not work with table substitution
        elements = default_table()

        # Fill in density or cell_volume.
        M = parse_formula(formula, natural_density=density)
        # CRUFT: use of T rather than H[1] is deprecated since 1.5.3
        if elements.T in M.atoms:
            warnings.warn("Use of tritium for labile hydrogen is deprecated."
                          " Use H[1] instead of T in your formula.")
            M = M.replace(elements.T, elements.H[1])
        if cell_volume is not None:
            # Note: cell_volume is only zero if there are no components
            M.density = 1e24*M.molecular_mass/cell_volume if cell_volume > 0 else 0
            #print name, M.molecular_mass, cell_volume, M.density
        else:
            cell_volume = 1e24*M.molecular_mass/M.density

        H = M.replace(elements.H[1], elements.H)
        D = M.replace(elements.H[1], elements.D)

        self.name = name
        self.cell_volume = cell_volume
        self.sld, self.Dsld = neutron_sld(H)[0], neutron_sld(D)[0]
        self.mass, self.Dmass = H.mass, D.mass
        self.D2Omatch = D2Omatch(self.sld, self.Dsld)
        self.charge = charge
        self.natural_formula = H
        self.labile_formula = M

        # TODO: formula should be natural_formula to be consistent
        # with sld and mass, which are computed with H-substitution.
        self.formula = self.labile_formula

    # TODO: are sld values float or complex?
    def D2Osld(self, volume_fraction: float=1., D2O_fraction: float=0.) -> float:
        """
        Neutron SLD of the molecule in a deuterated solvent.

        Changed 1.5.3: fix errors in SLD calculations.
        """
        solvent_sld = D2O_fraction*D2O_SLD + (1-D2O_fraction)*H2O_SLD
        solute_sld = D2O_fraction*self.Dsld + (1-D2O_fraction)*self.sld
        return volume_fraction*solute_sld + (1-volume_fraction)*solvent_sld

class Sequence(Molecule):
    """
    Convert FASTA sequence into chemical formula.

    *name* sequence name

    *sequence* code string

    *type* is one of::

       aa: amino acid sequence
       dna: dna sequence
       rna: rna sequence

    Note: rna sequence files treat T as U and dna sequence files treat U as T.
    """
    sequence: str

    @staticmethod
    def loadall(filename: Path|str, type: str|None=None) -> Iterator["Sequence"]:
        """
        Iterate over sequences in FASTA file, loading each in turn.

        Yields one FASTA sequence each cycle.
        """
        type = _guess_type_from_filename(str(filename), type)
        with open(filename, 'rt') as fh:
            for name, seq in read_fasta(fh):
                yield Sequence(name, seq, type=type)

    @staticmethod
    def load(filename: Path|str, type=None) -> "Sequence":
        """
        Load the first FASTA sequence from a file.
        """
        type = _guess_type_from_filename(str(filename), type)
        with open(filename, 'rt') as fh:
            name, seq = next(read_fasta(fh))
            return Sequence(name, seq, type=type)

    def __init__(self, name: str, sequence: str, type: str='aa'):
        # TODO: duplicated in Molecule.__init__
        # TODO: fasta does not work with table substitution
        elements = default_table()

        codes = CODE_TABLES[type]
        sequence = sequence.split('*', 1)[0]  # stop at first '*'
        sequence = sequence.replace(' ', '')  # ignore spaces
        parts = tuple(codes[c] for c in sequence)
        cell_volume = sum(p.cell_volume for p in parts)
        charge = sum(p.charge for p in parts)
        structure = []
        for p in parts:
            structure.extend(list(p.labile_formula.structure))
        # Add H + OH terminators to the sequence
        structure.extend(((2, elements.H[1]), (1, elements.O)))
        formula = parse_formula(structure).hill

        Molecule.__init__(
            self, name, formula, cell_volume=cell_volume, charge=charge)
        self.sequence = sequence

def _guess_type_from_filename(filename: str, type: str|None) -> str:
    if type is None:
        if filename.endswith('.fna'):
            type = 'dna'
        elif filename.endswith('.ffn'):
            type = 'dna'
        elif filename.endswith('.faa'):
            type = 'aa'
        elif filename.endswith('.frn'):
            type = 'rna'
        else:
            type = 'aa'
    return type

# PAK: Fixed in 1.5.3. Previous calculation used H2O density rather than D2O.
#: real portion of H2O sld at 20 C
H2O_SLD = neutron_sld("H2O@0.9982n")[0]
#: real portion of D2O sld at 20 C
#: Change 1.5.2: Use correct density in SLD calculation
D2O_SLD = neutron_sld("D2O@0.9982n")[0]
def D2Omatch(Hsld: float, Dsld: float) -> float:
    """
    Find the D2O% concentration of solvent such that neutron SLD of the
    material matches the neutron SLD of the solvent.

    *Hsld*, *Dsld* are the SLDs for the hydrogenated and deuterated forms
    of the material respectively, where *D* includes all the labile protons
    swapped for deuterons.  Water SLD is calculated at 20 C.

    Note that the resulting percentage is only meaningful between
    0% to 100%.  Beyond 100% you will need an additional constrast agent
    in the 100% D2O solvent to increase the SLD enough to match.

    .. deprecated:: 1.5.3
        Use periodictable.nsf.D2O_match(formula) instead.

    Change 1.5.3: corrected D2O sld, which will change the computed match point.
    """
    # SLD(%Dsample + (1-%)Hsample) = SLD(%D2O + (1-%)H2O)
    # => %SLD(Dsample) + (1-%)SLD(Hsample) = %SLD(D2O) + (1-%)SLD(H2O)
    # => %(SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O))
    #       = SLD(H2O) - SLD(Hsample)
    # => % = 100*(SLD(H2O) - SLD(Hsample))
    #           / (SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O))
    return 100 * (H2O_SLD - Hsld) / (Dsld - Hsld + H2O_SLD - D2O_SLD)


def read_fasta(fp: IO[str]) -> Iterator[tuple[str, str]]:
    """
    Iterate over the sequences in a FASTA file.

    Each iteration is a pair (sequence name, sequence codes).

    Change 1.5.3: Now uses H[1] rather than T for labile hydrogen.
    """
    name = ""
    seq: list[str] = []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield name, ''.join(seq)

def _code_average(bases, code_table) -> tuple[Formula, float, float]:
    """
    Compute average over possible nucleotides, assuming equal weight if
    precise nucleotide is not known.

    Note: averaging can lead to a fractional charge on the returned molecule.
    """
    n = len(bases)
    formula, cell_volume, charge = parse_formula(), 0., 0.
    for c in bases:
        base = code_table[c]
        formula += base.labile_formula
        cell_volume += base.cell_volume
        charge += base.charge
    if n > 0:
        formula, cell_volume, charge = (1/n) * formula, cell_volume/n, charge/n
    return formula, cell_volume, charge

def _set_amino_acid_average(target: str, codes: str, name: str|None=None) -> None:
    """
    Fill in partial unknowns for amino acids, such as "B" for aspartic acid or asparagine.
    """
    formula, cell_volume, charge = _code_average(codes, AMINO_ACID_CODES)
    if name is None:
        name = "/".join(AMINO_ACID_CODES[c].name for c in codes)
    molecule = Molecule(name, formula, cell_volume=cell_volume, charge=charge)
    AMINO_ACID_CODES[target] = molecule

# TODO: importing fasta does work, computing the neutron SLD for each molecule.
# This triggers nsf.init() which defines the neutron data for each isotope.
# Further, this does not allow private tables for fasta calculations.

# FASTA code table
def _(code: str, V: float, formula: str, name: str) -> tuple[str, Molecule]:
    if formula[-1] == '-':
        charge = -1
        formula = formula[:-1]
    elif formula[-1] == '+':
        charge = +1
        formula = formula[:-1]
    else:
        charge = 0
    molecule = Molecule(name, formula, cell_volume=V, charge=charge)
    molecule.code = code  # Add code attribute so we can write as well as read
    return code, molecule

# pylint: disable=bad-whitespace
AMINO_ACID_CODES: dict[str, Molecule] = dict((
    #code, volume, formula,        name
    _("A",  91.5, "C3H4H[1]NO",    "alanine"),
    #B: D or N
    _("C", 105.6, "C3H3H[1]NOS",   "cysteine"),
    _("D", 124.5, "C4H3H[1]NO3-",  "aspartic acid"),
    _("E", 155.1, "C5H5H[1]NO3-",  "glutamic acid"),
    _("F", 203.4, "C9H8H[1]NO",    "phenylalanine"),
    _("G",  66.4, "C2H2H[1]NO",    "glycine"),
    _("H", 167.3, "C6H5H[1]3N3O+", "histidine"),
    _("I", 168.8, "C6H10H[1]NO",   "isoleucine"),
    #J: L or I
    _("K", 171.3, "C6H9H[1]4N2O+", "lysine"),
    _("L", 168.8, "C6H10H[1]NO",   "leucine"),
    _("M", 170.8, "C5H8H[1]NOS",   "methionine"),
    _("N", 135.2, "C4H3H[1]3N2O2", "asparagine"),
    #O: _("O", ???.?, "C12H21N3O3", "pyrrolysine") -- update X below
    _("P", 129.3, "C5H7NO",     "proline"),
    _("Q", 161.1, "C5H5H[1]3N2O2", "glutamine"),
    _("R", 202.1, "C6H7H[1]6N4O+", "arginine"),
    _("S",  99.1, "C3H3H[1]2NO2",  "serine"),
    _("T", 122.1, "C4H5H[1]2NO2",  "threonine"),
    #U: selenocysteine -- update X below
    _("V", 141.7, "C5H8H[1]NO",    "valine"),
    _("W", 237.6, "C11H8H[1]2N2O", "tryptophan"),
    #X: any
    _("Y", 203.6, "C9H7H[1]2NO2",  "tyrosine"),
    #Z: E or Q
    #-: gap
    ))
_set_amino_acid_average('B', 'DN')
_set_amino_acid_average('J', 'LI')
_set_amino_acid_average('Z', 'EQ')
_set_amino_acid_average('X', 'ACDEFGHIKLMNPQRSTVWY', name='any')
_set_amino_acid_average('-', '', name='gap')
__doc__ += "\n\n*AMINO_ACID_CODES*::\n\n    " + "\n    ".join(
    "%s: %s"%(k, v.name) for k, v in sorted(AMINO_ACID_CODES.items()))

# mypy doesn't like redefinitions
def _1(formula: str, V: float, name: str) -> tuple[str, Molecule]:
    molecule = Molecule(name, formula, cell_volume=V)
    return name, molecule
NUCLEIC_ACID_COMPONENTS: dict[str, Molecule] = dict((
    # formula, volume, name
    _1("NaPO3",      60, "phosphate"),
    _1("C5H6H[1]O3",   125, "ribose"),
    _1("C5H7O2",    115, "deoxyribose"),
    _1("C5H2H[1]2N5",  114, "adenine"),
    _1("C4H2H[1]N2O2",  99, "uracil"),
    _1("C5H4H[1]N2O2", 126, "thymine"),
    _1("C5HH[1]3N5O",  119, "guanine"),
    _1("C4H2H[1]2N3O", 103, "cytosine"),
    ))
__doc__ += "\n\n*NUCLEIC_ACID_COMPONENTS*::\n\n  " + "\n  ".join(
    "%s: %s"%(k, v.formula) for k, v in sorted(NUCLEIC_ACID_COMPONENTS.items()))

CARBOHYDRATE_RESIDUES: dict[str, Molecule] = dict((
    # formula, volume, name
    _1("C6H7H[1]3O5",    171.9, "Glc"),
    _1("C6H7H[1]3O5",    166.8, "Gal"),
    _1("C6H7H[1]3O5",    170.8, "Man"),
    _1("C6H7H[1]4O5",    170.8, "Man (terminal)"),
    _1("C8H10H[1]3NO5",  222.0, "GlcNAc"),
    _1("C8H10H[1]3NO5",  232.9, "GalNAc"),
    _1("C6H7H[1]3O4",    160.8, "Fuc (terminal)"),
    _1("C11H11H[1]5NO8", 326.3, "NeuNac (terminal)"),
    # Glycosaminoglycans
    _1("C14H15H[1]5NO11Na", 390.7, "hyaluronate"),  # GlcA.GlcNAc
    _1("C14H17H[1]5NO13SNa", 473.5, "keratan sulphate"), # Gal.GlcNAc.SO4
    _1("C14H15H[1]4NO14SNa", 443.5, "chondroitin sulphate"), # GlcA.GalNAc.SO4
    ))
__doc__ += "\n\n*CARBOHYDRATE_RESIDUES*::\n\n  " + "\n  ".join(
    "%s: %s"%(k, v.formula) for k, v in sorted(CARBOHYDRATE_RESIDUES.items()))

LIPIDS: dict[str, Molecule] = dict((
    # formula, volume, name
    _1("CH2", 27, "methylene"),
    _1("CD2", 27, "methylene-D"),
    _1("C10H18NO8P", 350, "phospholipid headgroup"),
    _1("C6H5O6", 240, "triglyceride headgroup"),
    _1("C36H72NO8P", 1089, "DMPC"),
    _1("C36H20D52NO8P", 1089, "DMPC-D52"),
    _1("C29H55H[1]3NO8P", 932, "DLPE"),
    _1("C27H45H[1]O", 636, "cholesteral"),
    _1("C45H78O2", 1168, "oleate"),
    _1("C57H104O6", 1617, "trioleate form"),
    _1("C39H77H[1]2N2O2P", 1166, "palmitate ester"),
    ))
__doc__ += "\n\n*LIPIDS*::\n\n  " + "\n  ".join(
    "%s: %s"%(k, v.formula) for k, v in sorted(LIPIDS.items()))


RNA_BASES: dict[str, Molecule] = {}
RNA_CODES: dict[str, Molecule] = {}
DNA_BASES: dict[str, Molecule] = {}
DNA_CODES: dict[str, Molecule] = {}

def _set_rna_dna_codes() -> None:
    """
    Convert RNA/DNA table values into Molecule.

    Measured volumes from isolated mers reported in Durchschlag 1997, converted
    from mL/mol, with an addition of 30.39 to account for phosphorylation and
    dehydration in sequence. The value of 30.39 comes from comparison of the
    volume given in Harroun (2006) to the volumes of the RNA + T nucleosides
    given in Buckin (1989) after correcting for units. Harroun (2006) doesn't
    give volumes for AGC in the DNA nucleosides despite them being different in
    the Buckin source (especially guanosine).
    """
    # code, formula, volume (mL/mol), name
    rna_bases = (
        ("A",  "C10H8H[1]3N5O6P", 170.8, "adenosine"),
        ("T",   "C9H8H[1]2N2O8P", 151.7, "uridine"), # Use H[1] for U in RNA
        ("G",  "C10H7H[1]4N5O7P", 178.2, "guanosine"),
        ("C",   "C9H8H[1]3N3O7P", 153.7, "cytidine"),
    )
    dna_bases = (
        ("A",  "C10H9H[1]2N5O5P", 169.8, "adenosine"),
        ("T", "C10H11H[1]1N2O7P", 167.6, "thymidine"),
        ("G",  "C10H8H[1]3N5O6P", 173.7, "guanosine"),
        ("C",   "C9H9H[1]2N3O6P", 153.4, "cytidine"),
    )

    codes = (
        #code, nucleotides,  name
        ("A", "A",    "adenosine"),
        ("C", "C",    "cytidine"),
        ("G", "G",    "guanosine"),
        ("T", "T",    "thymidine"),
        ("U", "T",    "uridine"), # RNA_BASES["T"] is uridine
        ("R", "AG",   "purine"),
        ("Y", "CT",   "pyrimidine"),
        ("K", "GT",   "ketone"),
        ("M", "AC",   "amino"),
        ("S", "CG",   "strong"),
        ("W", "AT",   "weak"),
        ("B", "CGT",  "not A"),
        ("D", "AGT",  "not C"),
        ("H", "ACT",  "not G"),
        ("V", "ACG",  "not T"),
        ("N", "ACGT", "any base"),
        ("X", "",     "masked"),
        ("-", "",     "gap"),
        )

    for code, formula, volume, name in rna_bases:
        cell_volume = volume * 1e24/avogadro_number + 30.39
        molecule = Molecule(name, formula, cell_volume=cell_volume)
        molecule.code = code
        RNA_BASES[code] = molecule

    for code, formula, volume, name in dna_bases:
        cell_volume = volume * 1e24/avogadro_number + 30.39
        molecule = Molecule(name, formula, cell_volume=cell_volume)
        molecule.code = code
        DNA_BASES[code] = molecule

    for code, bases, name in codes:
        D, V, _ = _code_average(bases, RNA_BASES)
        rna = Molecule(name, D.hill, cell_volume=V)
        rna.code = code
        D, V, _ = _code_average(bases, DNA_BASES)
        dna = Molecule(name, D.hill, cell_volume=V)
        dna.code = code
        RNA_CODES[code] = rna
        DNA_CODES[code] = dna

_set_rna_dna_codes()
# pylint: enable=bad-whitespace

__doc__ += "\n\n*RNA_BASES*::\n\n  " + "\n  ".join(
    "%s:%s"%(k, v.name) for k, v in sorted(RNA_BASES.items()))
__doc__ += "\n\n*DNA_BASES*::\n\n  " + "\n  ".join(
    "%s:%s"%(k, v.name) for k, v in sorted(DNA_BASES.items()))


CODE_TABLES: dict[str, dict[str, Molecule]] = {
    'aa': AMINO_ACID_CODES,
    'dna': DNA_CODES,
    'rna': RNA_CODES,
}

def fasta_table() -> None:
    elements = default_table()

    rows = []
    rows += [v for k, v in sorted(AMINO_ACID_CODES.items())]
    rows += [v for k, v in sorted(NUCLEIC_ACID_COMPONENTS.items())]
    rows += [Sequence("beta casein", beta_casein)]

    print("%25s %7s %7s %7s %5s %5s %5s %5s %5s %5s"
          % ("name", "M(H2O)", "M(D2O)", "volume",
             "den", "#el", "xray", "nH2O", "nD2O", "%D2O match"))
    for v in rows:
        protons = sum(num*el.number for el, num in v.natural_formula.atoms.items())
        electrons = protons - v.charge
        Xsld = xray_sld(v.formula, wavelength=elements.Cu.K_alpha)
        print("%25s %7.1f %7.1f %7.1f %5.2f %5d %5.2f %5.2f %5.2f %5.1f"%(
            v.name, v.mass, v.Dmass, v.cell_volume, v.natural_formula.density,
            electrons, Xsld[0], v.sld, v.Dsld, v.D2Omatch))

beta_casein = "RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV"

## Uncomment to show package path on CI infrastructure
#def doctestpath():
#    """
#    Checking import path for doctests::
#
#        >>> import periodictable
#        >>> print(f"Path to imported periodictable in docstr is {periodictable.__file__}")
#        some path printed here
#    """

def test():
    from periodictable.constants import avogadro_number
    from .formulas import formula
    elements = default_table()

    ## Uncomment to show package path on CI infrastructure
    #import periodictable
    #print(f"Path to imported periodictable in package is {periodictable.__file__}")
    #print(fail_test)

    # Beta casein results checked against Duncan McGillivray's spreadsheet
    # name        Hmass   Dmass   vol     den   #el   xray  Hsld  Dsld
    # =========== ======= ======= ======= ===== ===== ===== ===== =====
    # beta casein 23561.9 23880.9 30872.9  1.27 12614 11.55  1.68  2.75
    # ... updated for new mass table [2023-08]
    #   same      23562.3 23881.2   same   1.27 same         1.68  2.75

    seq = Sequence("beta casein", beta_casein)
    density = seq.mass/avogadro_number/seq.cell_volume*1e24
    # print(seq.formula)
    # print(seq.mass, seq.Dmass, density, seq.sld, seq.Dsld)

    # Sequence now includes terminators: H[1]-...-OH[1], so adjust the masses
    # slightly from the spreadsheet values.
    H2O = formula("H2O@1n")
    D2O = formula("D2O@1n")
    assert abs(seq.mass - (23562.3 + H2O.mass)) < 0.1
    assert abs(seq.Dmass - (23881.2 + D2O.mass)) < 0.1
    assert abs(seq.cell_volume - 30872.9) < 0.1
    assert abs(density - 1.267) < 0.01
    assert abs(seq.sld - 1.68) < 0.01
    assert abs(seq.Dsld - 2.75) < 0.01

    # Check that X-ray sld is independent of isotope
    H = seq.labile_formula.replace(elements.H[1], elements.H)
    D = seq.labile_formula.replace(elements.H[1], elements.D)
    Hsld, Dsld = xray_sld(H, wavelength=1.54), xray_sld(D, wavelength=1.54)
    #print Hsld, Dsld
    assert abs(Hsld[0]-Dsld[0]) < 1e-10

if __name__ == "__main__":
    fasta_table()
    #test()
