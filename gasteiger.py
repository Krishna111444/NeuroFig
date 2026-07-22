"""
gasteiger.py — Gasteiger-Marsili partial atomic charges from first principles.

This is a from-scratch implementation of the PEOE method (Partial Equalization of
Orbital Electronegativity; Gasteiger & Marsili, Tetrahedron 1980, 36, 3219) — the
exact quantity that was silently zeroed in our PDBQT bug.

The science, in three lines:
  • each atom's electronegativity depends on its own charge:  χ(Q) = a + b·Q + c·Q²
  • across every bond, charge flows toward the more electronegative atom, an amount
    (χ_acceptor − χ_donor) / χ⁺_donor, where χ⁺ is the donor's cation electronegativity
  • damp the transfer by (1/2)^k each iteration and repeat ~6× until it converges

It is a deterministic iterative algorithm — no machine learning, no training. The
element/hybridisation parameters below are the published Gasteiger-Marsili values;
correctness is checked by reproducing Open Babel's charges (see the validation
script and tests/test_gasteiger.py).
"""
from __future__ import annotations

# (a, b, c) in eV, per element and hybridisation state. χ⁺ = a + b + c (cation
# electronegativity), except hydrogen which uses the special value 20.02.
PARAMS: dict[str, tuple[float, float, float]] = {
    "H":   (7.17, 6.24, -0.56),
    "C.3": (7.98, 9.18, 1.88),   # sp3 carbon
    "C.2": (8.79, 9.32, 1.51),   # sp2 / aromatic carbon
    "C.1": (10.39, 9.45, 0.73),  # sp carbon
    "N.3": (11.54, 10.82, 1.36),
    "N.2": (12.87, 11.15, 0.85),
    "N.1": (15.68, 11.70, -0.27),
    "O.3": (14.18, 12.92, 1.39),
    "O.2": (17.07, 13.79, 0.47),
    "F":   (14.66, 13.85, 2.31),
    "Cl":  (11.00, 9.69, 1.35),
    "Br":  (10.08, 8.47, 1.16),
    "I":   (9.90, 7.96, 0.96),
    "S.3": (10.14, 9.13, 1.38),
    "P.3": (8.90, 8.32, 1.48),
}
_CHI_PLUS_H = 20.02  # hydrogen's cation electronegativity is fixed high by the method


def chi_plus(param_key: str) -> float:
    if param_key == "H":
        return _CHI_PLUS_H
    a, b, c = PARAMS[param_key]
    return a + b + c


def sybyl_to_key(sybyl: str) -> str:
    """Map a SYBYL atom type (from a MOL2 file) to a parameter key."""
    el = sybyl.split(".")[0]
    if el == "H":
        return "H"
    if el == "C":
        return "C.1" if sybyl == "C.1" else ("C.3" if sybyl == "C.3" else "C.2")
    if el == "N":
        return "N.1" if sybyl == "N.1" else ("N.3" if sybyl in ("N.3", "N.4") else "N.2")
    if el == "O":
        return "O.3" if sybyl == "O.3" else "O.2"
    if el in ("F", "Cl", "Br", "I"):
        return el
    if el == "S":
        return "S.3"
    if el == "P":
        return "P.3"
    raise KeyError(f"no Gasteiger parameters for atom type {sybyl!r}")


def _element(key: str) -> str:
    return key.split(".")[0]


def resonance_retype(param_keys: list[str], bonds: list[tuple[int, int]],
                     bond_orders: list) -> list[str]:
    """Give conjugated heteroatoms sp2 character before charge equalisation.

    In a carboxyl/ester/amide, the single-bonded O or N is not really sp3 — its
    lone pair delocalises into the adjacent C=O (resonance), so it behaves as sp2.
    A plain PEOE pass misses this and under-charges the atom by ~0.15 e; retyping
    such atoms to sp2 reproduces Open Babel exactly. Rule: an sp3 O/N joined by a
    single bond to a carbon that also has a double bond to O/N becomes sp2.
    """
    keys = list(param_keys)
    carbonyl_c = set()
    for (a, b), t in zip(bonds, bond_orders):
        if str(t) == "2":
            for x, y in ((a, b), (b, a)):
                if _element(keys[x]) == "C" and _element(keys[y]) in ("O", "N"):
                    carbonyl_c.add(x)
    for (a, b), t in zip(bonds, bond_orders):
        if str(t) == "1":
            for x, y in ((a, b), (b, a)):
                if keys[x] == "O.3" and y in carbonyl_c:
                    keys[x] = "O.2"
                elif keys[x] == "N.3" and y in carbonyl_c:
                    keys[x] = "N.2"
    return keys


def compute_charges(param_keys: list[str], bonds: list[tuple[int, int]],
                    bond_orders: list | None = None, resonance: bool = True,
                    n_iter: int = 6) -> list[float]:
    """PEOE charges for atoms described by `param_keys` and a `bonds` edge list.

    param_keys  : parameter key per atom (e.g. 'C.3', 'H', 'O.2'), same order as atoms.
    bonds       : (i, j) index pairs (0-based).
    bond_orders : optional per-bond order ('1'/'2'/'ar'/2 …). When given (and
                  resonance=True), conjugated O/N are retyped to sp2 first — needed
                  for correct carboxyl/ester/amide charges. PEOE itself ignores order.
    """
    if bond_orders is not None and resonance:
        param_keys = resonance_retype(param_keys, bonds, bond_orders)
    n = len(param_keys)
    q = [0.0] * n
    p = [PARAMS[k] for k in param_keys]
    cp = [chi_plus(k) for k in param_keys]

    def chi(i):
        a, b, c = p[i]
        return a + b * q[i] + c * q[i] * q[i]

    for k in range(1, n_iter + 1):
        scale = 0.5 ** k
        dq = [0.0] * n
        for i, j in bonds:
            ci, cj = chi(i), chi(j)
            if ci >= cj:            # i is the acceptor, j the donor
                d = (ci - cj) / cp[j] * scale
                dq[i] -= d
                dq[j] += d
            else:                    # j is the acceptor, i the donor
                d = (cj - ci) / cp[i] * scale
                dq[j] -= d
                dq[i] += d
        for a in range(n):
            q[a] += dq[a]
    return q
