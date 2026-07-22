"""
Tests for gasteiger — the from-scratch PEOE charge calculator.

Two kinds of check:
  • backend-free regression values (methane, water) — the algorithm is
    deterministic, so exact numbers are pinned.
  • cross-validation against Open Babel when it is installed: molecules WITHOUT
    resonance edge cases must match OB to <0.005 e per atom. (Carboxylic acids
    differ because OB applies special resonance handling — documented, not hidden.)
Run: python -m tests.test_gasteiger
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gasteiger as g  # noqa: E402


def test_methane_regression():
    # C + 4 H, all C-H bonds. Carbon is slightly more electronegative than H,
    # so H go positive and C negative; total charge is ~0.
    keys = ["C.3", "H", "H", "H", "H"]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    q = g.compute_charges(keys, bonds)
    assert abs(sum(q)) < 1e-9                 # neutral molecule
    assert q[0] < 0 and all(h > 0 for h in q[1:])
    assert abs(q[0] - (-0.078)) < 0.002       # pinned value


def test_water_regression():
    keys = ["O.3", "H", "H"]
    bonds = [(0, 1), (0, 2)]
    q = g.compute_charges(keys, bonds)
    assert abs(sum(q)) < 1e-9
    assert abs(q[0] - (-0.410)) < 0.003       # oxygen strongly negative
    assert abs(q[1] - 0.205) < 0.003


def test_convergence_is_stable():
    keys = ["C.3", "H", "H", "H", "H"]
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4)]
    q6 = g.compute_charges(keys, bonds, n_iter=6)
    q12 = g.compute_charges(keys, bonds, n_iter=12)
    assert max(abs(a - b) for a, b in zip(q6, q12)) < 1e-3   # converged by 6 iters


def test_matches_openbabel_on_nonresonant():
    # Only runs if Open Babel is installed. Ethanol has no resonance edge case.
    import pdbqt_tools as pt
    if not pt.obabel_available():
        return
    import subprocess
    import tempfile
    exe = pt.find_obabel()
    with tempfile.TemporaryDirectory() as d:
        m = os.path.join(d, "m.mol2")
        subprocess.run([exe, "-:CCO", "-O", m, "--gen3d", "-h", "--partialcharge", "gasteiger"],
                       capture_output=True, text=True)
        atoms, bonds, obq, sec = [], [], [], None
        for line in open(m):
            if line.startswith("@<TRIPOS>"):
                sec = line.strip(); continue
            if sec == "@<TRIPOS>ATOM" and line.strip():
                f = line.split(); atoms.append(f[5]); obq.append(float(f[8]))
            elif sec == "@<TRIPOS>BOND" and line.strip():
                f = line.split(); bonds.append((int(f[1]) - 1, int(f[2]) - 1))
    mine = g.compute_charges([g.sybyl_to_key(a) for a in atoms], bonds)
    assert max(abs(mine[i] - obq[i]) for i in range(len(atoms))) < 0.005


def _mol2_from_smiles(smi):
    import pdbqt_tools as pt
    import subprocess
    import tempfile
    exe = pt.find_obabel()
    with tempfile.TemporaryDirectory() as d:
        m = os.path.join(d, "m.mol2")
        subprocess.run([exe, "-:" + smi, "-O", m, "--gen3d", "-h", "--partialcharge", "gasteiger"],
                       capture_output=True, text=True)
        atoms, bonds, orders, obq, sec = [], [], [], [], None
        for line in open(m):
            if line.startswith("@<TRIPOS>"):
                sec = line.strip(); continue
            if sec == "@<TRIPOS>ATOM" and line.strip():
                f = line.split(); atoms.append(f[5]); obq.append(float(f[8]))
            elif sec == "@<TRIPOS>BOND" and line.strip():
                f = line.split(); bonds.append((int(f[1]) - 1, int(f[2]) - 1)); orders.append(f[3])
    return atoms, bonds, orders, obq


def test_resonance_matches_openbabel_on_carboxyls():
    # Carboxylic acids / esters: resonance handling must reproduce OB exactly,
    # where a plain PEOE pass is off by ~0.15 e at the -OH oxygen.
    import pdbqt_tools as pt
    if not pt.obabel_available():
        return
    for smi in ("CC(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(=O)OC"):  # acetic, ibuprofen, ester
        atoms, bonds, orders, obq = _mol2_from_smiles(smi)
        keys = [g.sybyl_to_key(a) for a in atoms]
        mine = g.compute_charges(keys, bonds, bond_orders=orders)
        assert max(abs(mine[i] - obq[i]) for i in range(len(atoms))) < 0.005, smi
        # and without resonance it is clearly worse — proving the fix matters
        plain = g.compute_charges(keys, bonds, resonance=False)
        assert max(abs(plain[i] - obq[i]) for i in range(len(atoms))) > 0.1, smi


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
