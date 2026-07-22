"""
Tests for labcalc — every calculator checked against a KNOWN reference value.
That known-answer discipline is what makes these safe to ship. Run:
    python -m tests.test_labcalc
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import labcalc as lc  # noqa: E402


def _close(a, b, tol=1e-2):
    return abs(a - b) <= tol * max(1.0, abs(b))


def test_solution_prep():
    assert _close(lc.mass_for_molarity(1, 1, 58.44), 58.44)      # 1 M NaCl, 1 L
    assert _close(lc.stock_volume_for_target(10, 1, 10), 1.0)    # 10x -> 1x in 10 mL
    try:
        lc.stock_volume_for_target(1, 5, 10)  # target > stock
        assert False
    except ValueError:
        pass


def test_nucleic_acids():
    assert _close(lc.dsdna_conc_ng_uL(1.0), 50.0)
    assert _close(lc.dna_ng_to_pmol(650, 1000), 1.0)
    assert _close(lc.dna_copies(1, 1000), 9.264e8, tol=2e-3)     # qPCR standard
    assert _close(lc.gc_content("GGCC"), 1.0)
    assert _close(lc.gc_content("ATGC"), 0.5)
    assert lc.reverse_complement("ATGC") == "GCAT"
    assert lc.translate("ATGGCCTGA") == "MA"                     # M, A, stop


def test_protein():
    assert _close(lc.protein_average_mw("G"), 75.07)             # glycine
    assert _close(lc.protein_average_mw("GG"), 132.12)           # glycylglycine
    assert lc.extinction_coeff_280("WY") == 6990                 # 5500 + 1490
    assert lc.extinction_coeff_280("WY", cystines=2) == 7240     # + 2*125


def test_spectro_and_buffer():
    assert _close(lc.beer_lambert_concentration(1.0, 50, 1), 0.02)
    assert _close(lc.henderson_hasselbalch_ph(4.76, 1, 1), 4.76)
    assert _close(lc.henderson_hasselbalch_ph(4.76, 10, 1), 5.76)


def test_micro_and_pharma():
    assert _close(lc.cfu_per_ml(50, 1e-6, 0.1), 5e8)
    assert _close(lc.doubling_time(1, 8, 6), 2.0)                # 3 doublings in 6 h
    assert _close(lc.human_equivalent_dose(10, "mouse"), 0.8108)  # FDA Km 3/37
    try:
        lc.human_equivalent_dose(10, "dragon")
        assert False
    except ValueError:
        pass


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
