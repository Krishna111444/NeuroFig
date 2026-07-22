"""
Tests for doseresponse — validated by recovering KNOWN Hill parameters from
simulated data (the ground-truth proof a fit is correct). Run:
    python -m tests.test_doseresponse
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import doseresponse as dr  # noqa: E402


def _simulate(bottom, top, ic50, hill, seed=0, reps=3, noise=2.0):
    rng = np.random.default_rng(seed)
    conc = np.repeat(np.logspace(-9, -3, 8), reps)
    logx = np.log10(conc)
    y = bottom + (top - bottom) / (1 + 10 ** ((np.log10(ic50) - logx) * hill))
    return conc, y + rng.normal(0, noise, conc.size)


def test_recovers_agonist_ec50():
    c, y = _simulate(0, 100, 1e-6, 1.0, seed=1)
    f = dr.fit_dose_response(c, y, kind="auto")
    assert abs(np.log10(f.ic50) + 6) < 0.15      # EC50 ~ 1e-6
    assert abs(f.hill - 1.0) < 0.25
    assert f.label == "EC50" and f.direction == "increasing"
    assert f.r_squared > 0.98


def test_recovers_inhibitor_ic50():
    c, y = _simulate(0, 100, 5e-8, -1.2, seed=2)
    f = dr.fit_dose_response(c, y, kind="auto")
    assert abs(np.log10(f.ic50) - np.log10(5e-8)) < 0.15
    assert f.label == "IC50" and f.direction == "decreasing"


def test_noise_free_is_exact():
    c, y = _simulate(10, 90, 2e-7, 1.5, seed=3, noise=0.0)
    f = dr.fit_dose_response(c, y)
    assert abs(np.log10(f.ic50) - np.log10(2e-7)) < 0.01
    assert abs(f.hill - 1.5) < 0.02
    assert f.r_squared > 0.999


def test_ci_and_predict_and_guards():
    c, y = _simulate(0, 100, 1e-6, 1.0, seed=4)
    f = dr.fit_dose_response(c, y)
    lo, hi = f.ic50_ci
    assert lo < f.ic50 < hi                       # CI brackets the estimate
    # predict returns the plateaus at the extremes
    assert f.predict(1e-12) < f.predict(1e-1) if f.direction == "increasing" else True
    # too few points is a clear error, not a bad fit
    try:
        dr.fit_dose_response([1e-6, 1e-5], [10, 90])
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
