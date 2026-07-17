"""
pdbqt_tools.py — prepare protein/ligand structures for AutoDock/Vina docking.

Two stages, deliberately separated:

1. **Clean the PDB** (pure Python, no dependencies): strip waters, keep chosen
   chains, drop alternate-location duplicate atoms, keep only the first model,
   optionally remove other heteroatoms, or pull a ligand out to its own file.
   These are the fiddly steps people get wrong by hand.

2. **Convert to PDBQT** (via Open Babel): add hydrogens, assign Gasteiger partial
   charges and AutoDock atom types. This *must* be done by a real chemistry
   engine — a hand-rolled converter would silently produce wrong charges/atom
   types and break docking. If Open Babel is not installed we say so clearly and
   never emit a bogus file.

Open Babel lives in its own conda env so it never fights the app's Python:
    conda create -n neurofig_chem -c conda-forge python=3.11 openbabel -y
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

# Residue names treated as solvent/water and removed by default.
WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "TIP3", "TIP4"}
_ATOM_RECORDS = ("ATOM  ", "HETATM")


@dataclass
class PDBSummary:
    n_atoms: int = 0
    n_hetatm: int = 0
    n_waters: int = 0
    chains: list[str] = field(default_factory=list)
    hetero_residues: list[str] = field(default_factory=list)  # non-water HETATM resnames
    n_models: int = 1


def _is_atom(line: str) -> bool:
    return line.startswith(_ATOM_RECORDS)


def _res_name(line: str) -> str:
    return line[17:20].strip().upper()


def _chain_id(line: str) -> str:
    return line[21] if len(line) > 21 else " "


def _alt_loc(line: str) -> str:
    return line[16] if len(line) > 16 else " "


def summarize_pdb(text: str) -> PDBSummary:
    """Report what a PDB contains so the user can choose what to keep."""
    s = PDBSummary(n_models=0)
    chains, heteros = [], []
    seen_first_model_done = False
    for line in text.splitlines():
        if line.startswith("MODEL"):
            s.n_models += 1
            continue
        if not _is_atom(line):
            continue
        # only describe the first model to avoid double counting NMR ensembles
        if s.n_models > 1 and seen_first_model_done:
            continue
        res = _res_name(line)
        ch = _chain_id(line)
        if ch not in chains and ch.strip():
            chains.append(ch)
        if line.startswith("HETATM"):
            s.n_hetatm += 1
            if res in WATER_RESIDUES:
                s.n_waters += 1
            elif res not in heteros:
                heteros.append(res)
        else:
            s.n_atoms += 1
    if s.n_models == 0:
        s.n_models = 1
    s.chains, s.hetero_residues = chains, heteros
    return s


def clean_pdb(text: str, *, keep_chains: list[str] | None = None,
              remove_waters: bool = True, remove_hetero: bool = False,
              keep_hetero: list[str] | None = None,
              drop_altloc: bool = True) -> tuple[str, dict]:
    """Return (cleaned_pdb_text, report).

    keep_chains  : keep only these chain IDs (None = all).
    remove_waters: drop HOH/WAT/... (default True).
    remove_hetero: drop all non-water HETATM (ligands, ions, cofactors)...
    keep_hetero  : ...except these residue names (e.g. keep a metal cofactor).
    drop_altloc  : keep only altLoc blank or 'A' (removes duplicate conformers).
    Only the first model of a multi-model (NMR) file is kept.
    """
    keep_hetero = {h.upper() for h in (keep_hetero or [])}
    kept: list[str] = []
    removed = {"waters": 0, "hetero": 0, "other_chains": 0, "altloc": 0}
    in_model, past_first_model = 0, False

    for line in text.splitlines():
        if line.startswith("MODEL"):
            in_model += 1
            if in_model > 1:
                past_first_model = True
            continue
        if line.startswith("ENDMDL"):
            continue
        if line.startswith("TER"):
            kept.append(line)
            continue
        if not _is_atom(line):
            continue  # drop HEADER/REMARK/CONECT/etc. -> minimal clean structure
        if past_first_model:
            continue

        res, ch, alt = _res_name(line), _chain_id(line), _alt_loc(line)
        is_water = res in WATER_RESIDUES
        is_het = line.startswith("HETATM")

        if remove_waters and is_water:
            removed["waters"] += 1; continue
        if keep_chains is not None and ch not in keep_chains:
            removed["other_chains"] += 1; continue
        if drop_altloc and alt not in (" ", "A"):
            removed["altloc"] += 1; continue
        if is_het and not is_water and remove_hetero and res not in keep_hetero:
            removed["hetero"] += 1; continue

        kept.append(line)

    if not kept or kept[-1].strip() != "END":
        kept.append("END")
    return "\n".join(kept) + "\n", removed


def extract_ligand(text: str, resname: str, chain: str | None = None) -> str:
    """Pull a single ligand (by residue name, optional chain) into its own PDB."""
    resname = resname.strip().upper()
    out = []
    for line in text.splitlines():
        if not _is_atom(line):
            continue
        if _res_name(line) == resname and (chain is None or _chain_id(line) == chain):
            out.append(line)
    if not out:
        raise ValueError(f"No atoms found for ligand residue {resname!r}"
                         + (f" in chain {chain}" if chain else ""))
    out.append("END")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------- Open Babel bridge
def find_obabel() -> str | None:
    """Locate the Open Babel executable: env var, PATH, then conda env guesses."""
    env = os.environ.get("NEUROFIG_OBABEL")
    if env and os.path.exists(env):
        return env
    on_path = shutil.which("obabel")
    if on_path:
        return on_path
    # common conda layouts (Windows puts CLI exes under Library/bin)
    roots = [sys.prefix, os.path.dirname(os.path.dirname(sys.prefix))]
    names = ["obabel.exe", "obabel"]
    subdirs = [os.path.join("envs", "neurofig_chem", "Library", "bin"),
               os.path.join("envs", "neurofig_chem", "bin"),
               os.path.join("Library", "bin"), "bin"]
    for root in roots:
        for sub in subdirs:
            for name in names:
                cand = os.path.join(root, sub, name)
                if os.path.exists(cand):
                    return cand
    return None


def obabel_available() -> bool:
    return find_obabel() is not None


def pdb_to_pdbqt(in_path: str, out_path: str, *, mode: str = "receptor",
                 ph: float | None = 7.4, obabel: str | None = None) -> str:
    """Convert a (cleaned) PDB to PDBQT with Open Babel.

    mode='receptor' -> rigid receptor (-xr), no rotatable bonds.
    mode='ligand'   -> flexible ligand (rotatable bonds detected by Open Babel).
    ph              -> add hydrogens for this pH (None = don't add H).
    Gasteiger partial charges and AutoDock atom types are assigned by Open Babel.
    Raises RuntimeError with install guidance if Open Babel is missing.
    """
    exe = obabel or find_obabel()
    if not exe:
        raise RuntimeError(
            "Open Babel not found. Install it in a dedicated env:\n"
            "  conda create -n neurofig_chem -c conda-forge python=3.11 openbabel -y\n"
            "or set NEUROFIG_OBABEL to the obabel executable path.")
    if mode not in ("receptor", "ligand"):
        raise ValueError("mode must be 'receptor' or 'ligand'.")

    cmd = [exe, in_path, "-O", out_path]
    if ph is not None:
        cmd += ["-p", str(ph)]
    if mode == "receptor":
        cmd += ["-xr"]  # rigid receptor
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.exists(out_path):
        raise RuntimeError(f"Open Babel failed (exit {proc.returncode}):\n{proc.stderr.strip()}")
    return out_path
