"""
Tests for standardcurve — validated by (a) recovering known parameters and
(b) back-calculation round-trips: a known concentration -> its signal -> back to
the concentration. Run: python -m tests.test_standardcurve
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import standardcurve as sc  # noqa: E402


def test_linear_recovery_and_roundtrip():
    conc = np.repeat([0, 25, 50, 100, 200, 400], 3).astype(float)
    sig = 0.002 * conc + 0.05
    f = sc.fit_standard_curve(conc, sig, "linear")
    assert abs(f.params["slope"] - 0.002) < 1e-6
    assert abs(f.params["intercept"] - 0.05) < 1e-6
    assert f.r_squared > 0.999
    # round-trip: predict then back-calculate returns the input
    assert abs(f.back_calculate(f.predict(123.0)) - 123.0) < 1e-6


def test_qpcr_efficiency():
    copies = np.repeat(np.logspace(1, 7, 7), 3)
    ct = -3.32 * np.log10(copies) + 40
    f = sc.fit_standard_curve(copies, ct, "qpcr")
    assert abs(f.params["slope"] - (-3.32)) < 1e-3
    assert 95 < f.efficiency * 100 < 105       # ~100% efficiency
    # back-calc a known point
    assert abs(np.log10(f.back_calculate(-3.32 * 3 + 40)) - 3) < 1e-6


def test_4pl_roundtrip_and_out_of_range():
    c = np.logspace(-3, 2, 10)
    y = 0.05 + (2.0 - 0.05) / (1 + 10 ** ((np.log10(1.0) - np.log10(c)) * 1.2))
    f = sc.fit_standard_curve(c, y, "4pl")
    assert f.r_squared > 0.999
    # known x=5 -> its signal -> back to ~5
    y5 = 0.05 + (2.0 - 0.05) / (1 + 10 ** ((np.log10(1.0) - np.log10(5.0)) * 1.2))
    assert abs(f.back_calculate(y5) - 5.0) < 0.05
    # a signal above the top plateau is not invertible -> NaN, never a fake value
    assert np.isnan(f.back_calculate(10.0))


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
