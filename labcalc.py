"""
labcalc.py — deterministic bench calculators for the life sciences.

Every function here is a small, exact calculation that researchers do constantly,
by hand or with a patchwork of unvalidated web tools, and get wrong under
deadline. They are pure and reference-checkable (see tests/test_labcalc.py), which
is exactly what makes them safe to sell: no AI, no guessing, just the right number.

Domains: solution prep, nucleic acids, protein, spectrophotometry, buffers,
microbiology, pharmacology. Grouped so the collection can grow one function at a
time without disturbing the rest.
"""
from __future__ import annotations

import math

AVOGADRO = 6.02214076e23
AVG_BP_MW = 650.0        # g/mol per base pair (double-stranded, common convention)
WATER_MW = 18.01528


# ---------------------------------------------------------------- solution prep
def mass_for_molarity(conc_mol_L: float, volume_L: float, mw_g_mol: float) -> float:
    """Grams of solute to make `volume_L` of a `conc_mol_L` solution. m = C·V·MW."""
    return conc_mol_L * volume_L * mw_g_mol


def stock_volume_for_target(stock_conc: float, target_conc: float,
                            target_volume: float) -> float:
    """Volume of stock needed for a dilution (C1·V1 = C2·V2) -> V1 = C2·V2/C1.

    Units are arbitrary but must match (conc with conc, volume with volume).
    """
    if stock_conc <= 0:
        raise ValueError("stock concentration must be > 0")
    if target_conc > stock_conc:
        raise ValueError("target concentration cannot exceed stock concentration")
    return target_conc * target_volume / stock_conc


# ---------------------------------------------------------------- nucleic acids
def dsdna_conc_ng_uL(a260: float, dilution_factor: float = 1.0) -> float:
    """dsDNA concentration from A260 (1.0 A260 = 50 ng/µL for dsDNA)."""
    return a260 * 50.0 * dilution_factor


def nucleic_purity(a260: float, a280: float) -> float:
    """A260/A280 purity ratio (~1.8 pure DNA, ~2.0 pure RNA)."""
    if a280 == 0:
        raise ValueError("A280 must be non-zero")
    return a260 / a280


def dna_ng_to_pmol(ng: float, length_bp: int) -> float:
    """Convert a mass of dsDNA (ng) to amount (pmol) given its length in bp."""
    if length_bp <= 0:
        raise ValueError("length_bp must be > 0")
    return ng * 1000.0 / (length_bp * AVG_BP_MW)


def dna_copies(ng: float, length_bp: int) -> float:
    """Molecule count in a mass of dsDNA — the qPCR standard-curve calculation."""
    if length_bp <= 0:
        raise ValueError("length_bp must be > 0")
    moles = (ng * 1e-9) / (length_bp * AVG_BP_MW)
    return moles * AVOGADRO


def gc_content(seq: str) -> float:
    """Fraction of G/C in a DNA sequence (0–1). Ignores non-ACGT characters."""
    s = [b for b in seq.upper() if b in "ACGT"]
    if not s:
        raise ValueError("no valid ACGT bases in sequence")
    return sum(b in "GC" for b in s) / len(s)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence (preserves case of the mapping)."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
            "a": "t", "t": "a", "g": "c", "c": "g", "n": "n"}
    return "".join(comp[b] for b in reversed(seq) if b in comp)


_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(seq: str, stop_at_stop: bool = True) -> str:
    """Translate a DNA coding sequence to one-letter protein (reads frame 1)."""
    s = seq.upper().replace(" ", "")
    out = []
    for i in range(0, len(s) - 2, 3):
        aa = _CODON_TABLE.get(s[i:i + 3], "X")
        if aa == "*" and stop_at_stop:
            break
        out.append(aa)
    return "".join(out)


# ---------------------------------------------------------------- protein
# Average residue masses (Da) — monomer mass minus one water (peptide bond).
_AA_AVG = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}


def protein_average_mw(seq: str) -> float:
    """Average molecular weight (Da) of a protein from its one-letter sequence."""
    s = [a for a in seq.upper() if a in _AA_AVG]
    if not s:
        raise ValueError("no valid amino acids in sequence")
    return sum(_AA_AVG[a] for a in s) + WATER_MW


def extinction_coeff_280(seq: str, cystines: int = 0) -> int:
    """Molar extinction coefficient at 280 nm (M⁻¹cm⁻¹), Gill & von Hippel.

    ε = nTrp·5500 + nTyr·1490 + nCystine·125. `cystines` = number of disulfide
    bonds (0 if all cysteines are reduced).
    """
    s = seq.upper()
    return s.count("W") * 5500 + s.count("Y") * 1490 + cystines * 125


# ---------------------------------------------------------------- spectrophotometry
def beer_lambert_concentration(absorbance: float, epsilon: float,
                               path_cm: float = 1.0) -> float:
    """Concentration from A = ε·c·l  ->  c = A/(ε·l). Units follow ε."""
    if epsilon <= 0 or path_cm <= 0:
        raise ValueError("epsilon and path length must be > 0")
    return absorbance / (epsilon * path_cm)


# ---------------------------------------------------------------- buffers
def henderson_hasselbalch_ph(pka: float, conc_base: float, conc_acid: float) -> float:
    """pH = pKa + log10([A⁻]/[HA]) for a weak-acid buffer."""
    if conc_base <= 0 or conc_acid <= 0:
        raise ValueError("concentrations must be > 0")
    return pka + math.log10(conc_base / conc_acid)


# ---------------------------------------------------------------- microbiology
def cfu_per_ml(colonies: int, dilution_factor: float, plated_ml: float) -> float:
    """Colony-forming units per mL = colonies / (dilution × volume plated)."""
    if dilution_factor <= 0 or plated_ml <= 0:
        raise ValueError("dilution and plated volume must be > 0")
    return colonies / (dilution_factor * plated_ml)


def doubling_time(n0: float, nt: float, elapsed: float) -> float:
    """Population doubling time = elapsed·ln2 / ln(Nt/N0) (units follow `elapsed`)."""
    if n0 <= 0 or nt <= 0 or nt <= n0 or elapsed <= 0:
        raise ValueError("need Nt > N0 > 0 and elapsed > 0")
    return elapsed * math.log(2) / math.log(nt / n0)


# ---------------------------------------------------------------- pharmacology
# FDA body-surface-area Km factors (mg/kg dose scaling across species).
KM_FACTORS = {"human": 37, "mouse": 3, "rat": 6, "hamster": 5, "guinea pig": 8,
              "rabbit": 12, "dog": 20, "monkey": 12, "cat": 10}


def human_equivalent_dose(animal_dose_mg_kg: float, species: str) -> float:
    """FDA human-equivalent dose (mg/kg) from an animal dose via BSA scaling.

    HED = animal_dose × (Km_animal / Km_human). Reference: FDA 2005 guidance.
    """
    sp = species.strip().lower()
    if sp not in KM_FACTORS:
        raise ValueError(f"unknown species {species!r}; known: {sorted(KM_FACTORS)}")
    return animal_dose_mg_kg * KM_FACTORS[sp] / KM_FACTORS["human"]
