"""
Inversion tests for the upload-driven figures. For each, we ask "how could this
break or mislead?" and feed exactly that input, asserting a scientific, non-
crashing result. Run: python -m tests.test_realdata
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import journal_figures as jf  # noqa: E402
import neurofig_core as nc  # noqa: E402


def _ok(fig):
    assert len(nc.figure_to_bytes(fig, "pdf")) > 1000
    plt.close(fig)


# ------------------------------------------------------------------ volcano
def test_volcano_sample_and_zero_p_floored():
    df = jf.sample_volcano_df()   # contains real p = 0 values
    fig, w = jf.volcano_from_data(df, "log2FoldChange", "pvalue", label_col="gene")
    _ok(fig)
    assert any("floored" in m for m in w)          # p=0 handled, not left as inf


def test_volcano_out_of_range_p_dropped():
    df = pd.DataFrame({"fc": [1, -2, 3, 0.5], "p": [0.01, 1.5, -0.2, 0.04]})
    fig, w = jf.volcano_from_data(df, "fc", "p")
    _ok(fig)
    assert any("outside [0, 1]" in m for m in w)    # 1.5 and -0.2 dropped


def test_volcano_linear_fc_nonpositive_dropped():
    df = pd.DataFrame({"fc": [2.0, 0.0, -1.0, 4.0], "p": [0.01, 0.02, 0.03, 0.04]})
    fig, w = jf.volcano_from_data(df, "fc", "p", fc_is_log2=False)
    _ok(fig)
    assert any("non-positive linear" in m for m in w)


def test_volcano_neglog10_mode_no_double_transform():
    # p already as -log10; must NOT be transformed again
    df = pd.DataFrame({"fc": [2, -2, 1], "nlp": [5.0, 4.0, 0.3]})
    fig, w = jf.volcano_from_data(df, "fc", "nlp", p_is_neglog10=True, p_thresh=2.0)
    _ok(fig)


def test_volcano_all_bad_raises():
    df = pd.DataFrame({"fc": [np.nan, np.nan], "p": [np.nan, np.nan]})
    try:
        jf.volcano_from_data(df, "fc", "p")
        assert False
    except ValueError:
        pass


# ------------------------------------------------------------------ raincloud
def test_raincloud_sample():
    fig, w = jf.raincloud_from_data(jf.sample_raincloud_df(), "group", "value")
    _ok(fig)


def test_raincloud_singleton_group_no_crash():
    df = pd.DataFrame({"g": ["A"] * 5 + ["B"], "v": [1.0, 2, 3, 4, 5, 9]})
    fig, w = jf.raincloud_from_data(df, "g", "v")   # B has n=1
    _ok(fig)
    assert any("points only" in m for m in w)


def test_raincloud_zero_variance_group_no_kde_crash():
    df = pd.DataFrame({"g": ["A"] * 6 + ["B"] * 6, "v": [7.0] * 6 + [1.0, 2, 3, 4, 5, 6]})
    fig, w = jf.raincloud_from_data(df, "g", "v")   # A is constant -> KDE singular
    _ok(fig)
    assert any("points only" in m for m in w)


def test_raincloud_nonnumeric_coerced():
    df = pd.DataFrame({"g": ["A", "A", "A", "B", "B", "B"],
                       "v": ["1", "2", "bad", "4", "5", "6"]})
    fig, w = jf.raincloud_from_data(df, "g", "v")
    _ok(fig)
    assert any("Dropped" in m for m in w)


# ------------------------------------------------------------------ correlation
def test_correlation_sample():
    fig, w = jf.correlation_heatmap_from_data(jf.sample_correlation_df(),
                                              list(jf.sample_correlation_df().columns))
    _ok(fig)


def test_correlation_constant_column_dropped():
    df = jf.sample_correlation_df().copy()
    df["flat"] = 3.0
    fig, w = jf.correlation_heatmap_from_data(df, list(df.columns))
    _ok(fig)
    assert any("constant column" in m for m in w)


def test_correlation_too_few_columns_raises():
    df = pd.DataFrame({"a": [1.0, 2, 3], "b": [5.0, 5, 5]})  # b constant -> only 1 left
    try:
        jf.correlation_heatmap_from_data(df, ["a", "b"])
        assert False
    except ValueError:
        pass


def test_correlation_handles_nan_and_spearman():
    df = jf.sample_correlation_df().copy()
    df.iloc[0, 0] = np.nan; df.iloc[5, 2] = np.nan
    fig, w = jf.correlation_heatmap_from_data(df, list(df.columns), method="spearman")
    _ok(fig)


def test_correlation_bad_method_raises():
    try:
        jf.correlation_heatmap_from_data(jf.sample_correlation_df(),
                                         list(jf.sample_correlation_df().columns),
                                         method="nonsense")
        assert False
    except ValueError:
        pass


def test_upset_needs_two_sets():
    df = pd.DataFrame({"A": ["g1", "g2"], "B": ["g2", "g3"]})
    fig, w = jf.upset_from_data(df, ["A", "B"])
    _ok(fig)
    try:
        jf.upset_from_data(df, ["A"]); assert False
    except ValueError:
        pass


def test_enrichment_floors_bad_padjust():
    df = pd.DataFrame({"term": ["T1", "T2", "T3"], "ratio": [0.1, 0.2, 0.3],
                       "count": [5, 10, 15], "padj": [0.0, 1e-4, 0.01]})  # p=0 present
    fig, w = jf.enrichment_from_data(df, "term", "ratio", "count", "padj")
    _ok(fig)
    assert any("floored" in m for m in w)


def test_pca_guards():
    import numpy as np
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.normal(0, 1, (20, 5)), columns=[f"f{i}" for i in range(5)])
    df["flat"] = 1.0  # constant -> dropped
    df["group"] = np.repeat(["a", "b"], 10)
    fig, w = jf.pca_from_data(df, [f"f{i}" for i in range(5)] + ["flat"], "group")
    _ok(fig)
    assert any("constant" in m for m in w)
    try:
        jf.pca_from_data(df, ["f0"]); assert False   # <2 features
    except ValueError:
        pass


def test_composition_drops_zero_sample_and_normalizes():
    df = pd.DataFrame({"sample": ["S1", "S2", "S3"],
                       "a": [1.0, 0.0, 2.0], "b": [3.0, 0.0, 2.0]})
    fig, w = jf.stacked_composition_from_data(df, "sample", ["a", "b"])
    _ok(fig)
    assert any("zero total" in m for m in w)   # S2 dropped


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
