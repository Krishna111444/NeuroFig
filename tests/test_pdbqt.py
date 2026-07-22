"""
Tests for pdbqt_tools — the pure-Python PDB-cleaning stage (fully validatable
without any chemistry backend) plus the missing-backend guard.

Run:  python -m tests.test_pdbqt   (or)  pytest -q
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pdbqt_tools as p  # noqa: E402

PDB = """HEADER    TEST
REMARK    dropped
ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.560  13.207  10.000  1.00  0.00           C
ATOM      3  CB AALA A   1      13.000  14.000  10.000  0.50  0.00           C
ATOM      4  CB BALA A   1      13.000  12.000  10.000  0.50  0.00           C
ATOM      5  N   GLY B   1      20.000  20.000  20.000  1.00  0.00           N
HETATM    6  O   HOH A 100      30.000  30.000  30.000  1.00  0.00           O
HETATM    7  C1  LIG A 200       5.000   5.000   5.000  1.00  0.00           C
HETATM    8  C2  LIG A 200       6.000   5.000   5.000  1.00  0.00           C
TER
END
"""


def _atom_lines(text):
    return [l for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]


def test_summary():
    s = p.summarize_pdb(PDB)
    assert s.n_atoms == 5 and s.n_hetatm == 3 and s.n_waters == 1
    assert s.chains == ["A", "B"] and s.hetero_residues == ["LIG"]


def test_clean_keep_chain_remove_water_hetero_altloc():
    clean, rep = p.clean_pdb(PDB, keep_chains=["A"], remove_waters=True,
                             remove_hetero=True, drop_altloc=True)
    lines = _atom_lines(clean)
    assert rep == {"waters": 1, "hetero": 2, "other_chains": 1, "altloc": 1}
    assert len(lines) == 3
    assert all(l[21] == "A" for l in lines)              # only chain A
    assert all(l[16] != "B" for l in lines)              # altLoc B removed
    assert all(l[17:20].strip() not in ("HOH", "LIG") for l in lines)


def test_keep_hetero_whitelist():
    clean, _ = p.clean_pdb(PDB, remove_hetero=True, keep_hetero=["LIG"])
    assert any(l[17:20].strip() == "LIG" for l in _atom_lines(clean))


def test_extract_ligand():
    lig = p.extract_ligand(PDB, "LIG")
    assert len(_atom_lines(lig)) == 2
    try:
        p.extract_ligand(PDB, "ZZZ")
        assert False, "expected ValueError for missing ligand"
    except ValueError:
        pass


def test_convert_guard_when_backend_missing():
    # If Open Babel is absent, we must raise a clear error, never emit junk.
    if p.obabel_available():
        return  # backend present: nothing to assert here
    try:
        p.pdb_to_pdbqt("in.pdb", "out.pdbqt")
        assert False, "expected RuntimeError when Open Babel missing"
    except RuntimeError as e:
        assert "Open Babel not found" in str(e)


def test_real_conversion_assigns_charges():
    # Only runs when Open Babel is installed. Verifies a real ligand converts to a
    # valid PDBQT WITH non-zero Gasteiger charges (zero charges break docking).
    if not p.obabel_available():
        return
    import subprocess
    import tempfile
    exe = p.find_obabel()
    with tempfile.TemporaryDirectory() as d:
        pdb = os.path.join(d, "lig.pdb"); out = os.path.join(d, "lig.pdbqt")
        subprocess.run([exe, "-:CC(C)Cc1ccc(cc1)C(C)C(=O)O", "-O", pdb, "--gen3d"],
                       capture_output=True, text=True)
        p.pdb_to_pdbqt(pdb, out, mode="ligand", ph=7.4)
        txt = open(out).read()
        assert "ROOT" in txt and ("BRANCH" in txt or "TORSDOF" in txt)
        atomlines = [l for l in txt.splitlines() if l.startswith(("ATOM", "HETATM"))]
        charges = [float(l[70:76]) for l in atomlines if len(l) >= 76 and l[70:76].strip()]
        assert any(abs(c) > 1e-6 for c in charges), "Gasteiger charges missing"


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
