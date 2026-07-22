"""
cluster_stats.py — cluster-based permutation test for comparing two time-courses.

The problem: comparing two conditions (e.g. CS+ vs CS-) across many time points
means many correlated tests. Testing each time bin separately does not control
the false-positive rate. The cluster-based permutation test (Maris & Oostenveld,
J. Neurosci. Methods 2007) solves this: it forms clusters of adjacent
supra-threshold time points, scores each by its summed statistic ("cluster mass"),
and calibrates significance against a permutation null of the *maximum* cluster
mass — which controls the family-wise error rate across time.

There is no clean standalone implementation for a plain (subjects x time) array
(FieldTrip/MNE tie it to their own data structures), so this is a real,
self-contained, reference-validatable component — not a wrapper.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def _clusters_from_tvals(t: np.ndarray, tcrit: float) -> list[tuple[int, int, float, int]]:
    """Contiguous runs of |t| > tcrit with a constant sign.

    Returns (start, end_inclusive, cluster_mass, sign) per cluster, where
    cluster_mass is the sum of t-values in the run.
    """
    supra = np.abs(t) > tcrit
    sign = np.sign(t)
    clusters, i, n = [], 0, len(t)
    while i < n:
        if supra[i]:
            j = i
            while j + 1 < n and supra[j + 1] and sign[j + 1] == sign[i]:
                j += 1
            clusters.append((i, j, float(np.sum(t[i:j + 1])), int(sign[i])))
            i = j + 1
        else:
            i += 1
    return clusters


def _max_abs_mass(clusters) -> float:
    return max((abs(m) for _, _, m, _ in clusters), default=0.0)


def cluster_permutation_test(A: np.ndarray, B: np.ndarray, *,
                             alpha: float = 0.05, cluster_alpha: float = 0.05,
                             n_perm: int = 1000, seed: int = 0) -> dict:
    """Independent-samples cluster permutation test between two conditions.

    A, B : arrays of shape (n_subjects, n_time). Same n_time; subject counts may
           differ. Rows are subjects/trials, columns are time points.
    cluster_alpha : pointwise threshold for forming clusters (two-sided t).
    alpha         : cluster-level significance.
    Returns dict with per-time t-values, the threshold, and the detected clusters,
    each carrying a permutation p-value and a significance flag.
    """
    A = np.asarray(A, float); B = np.asarray(B, float)
    nA, nt = A.shape
    nB = B.shape[0]
    df = nA + nB - 2
    tcrit = float(stats.t.ppf(1 - cluster_alpha / 2, df))

    t_obs = stats.ttest_ind(A, B, axis=0).statistic
    obs = _clusters_from_tvals(t_obs, tcrit)

    pooled = np.concatenate([A, B], axis=0)
    rng = np.random.default_rng(seed)
    null_max = np.empty(n_perm)
    for k in range(n_perm):
        idx = rng.permutation(nA + nB)
        tp = stats.ttest_ind(pooled[idx[:nA]], pooled[idx[nA:]], axis=0).statistic
        null_max[k] = _max_abs_mass(_clusters_from_tvals(tp, tcrit))

    clusters = []
    for (s, e, mass, sign) in obs:
        # unbiased permutation p-value: (1 + #{null >= observed}) / (n_perm + 1)
        p = (1 + int(np.sum(null_max >= abs(mass)))) / (n_perm + 1)
        clusters.append({"start": int(s), "end": int(e), "mass": mass,
                         "sign": sign, "p": float(p), "sig": bool(p < alpha)})

    return {"t": t_obs, "tcrit": tcrit, "clusters": clusters, "null_max": null_max}
