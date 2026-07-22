"""
Tests for cluster_stats — the cluster-based permutation test.

Beyond ordinary unit tests, we assert the *defining* property of a valid test:
false-positive control under the null. Run: python -m tests.test_cluster
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cluster_stats as cs  # noqa: E402


def test_cluster_finder_exact():
    t = np.array([0, 0, 3.0, 4, 3, 0, -3, -4, 0])
    assert cs._clusters_from_tvals(t, 2.0) == [(2, 4, 10.0, 1), (6, 7, -7.0, -1)]


def test_power_detects_and_localizes():
    rng = np.random.default_rng(1)
    nt = 100
    A = rng.normal(0, 1, (15, nt)); B = rng.normal(0, 1, (15, nt))
    A[:, 40:60] += 1.5  # true effect only here
    r = cs.cluster_permutation_test(A, B, n_perm=300, seed=7)
    sig = [c for c in r["clusters"] if c["sig"]]
    assert sig, "should detect a real effect"
    big = max(sig, key=lambda c: abs(c["mass"]))
    assert not (big["end"] < 40 or big["start"] >= 60), "cluster should overlap true window"


def test_false_positive_rate_controlled():
    # Under the null (same distribution), a valid test flags <= ~alpha of runs.
    # Small MC for speed; a broken test would blow well past the loose bound.
    n_sim, hits = 80, 0
    for s in range(n_sim):
        rng = np.random.default_rng(5000 + s)
        A = np.cumsum(rng.normal(0, 1, (10, 50)), axis=1) / np.sqrt(50)
        B = np.cumsum(rng.normal(0, 1, (10, 50)), axis=1) / np.sqrt(50)
        r = cs.cluster_permutation_test(A, B, n_perm=100, seed=s)
        hits += any(c["sig"] for c in r["clusters"])
    fpr = hits / n_sim
    assert fpr <= 0.15, f"false-positive rate {fpr:.2f} too high — test invalid"


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
