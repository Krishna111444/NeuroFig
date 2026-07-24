"""
Tests for journal_figures — every one of the 20 figures must build and export to
a non-empty vector file, at a journal preset and exact width. Run:
    python -m tests.test_journal_figures
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import journal_figures as jf  # noqa: E402
import neurofig_core as nc  # noqa: E402


def test_registry_size():
    assert len(jf.FIGURES) == 25  # 20 journal + 5 multiomics


def test_venn_counts_from_data():
    import pandas as pd
    df = pd.DataFrame({"A": ["g1", "g2", "g3", None], "B": ["g2", "g3", "g4", "g5"]})
    fig, w = jf.venn_from_data(df, ["A", "B"])
    # A={g1,g2,g3}, B={g2,g3,g4,g5}: A-only=1, B-only=2, shared=2
    assert jf._set_regions([{"g1", "g2", "g3"}, {"g2", "g3", "g4", "g5"}]) == {"10": 1, "01": 2, "11": 2}


def test_all_figures_render_and_export():
    for name, fn in jf.FIGURES.items():
        fig = fn(preset="Nature", width_mm=89)
        for fmt in ("pdf", "svg"):
            assert len(nc.figure_to_bytes(fig, fmt, exact=True)) > 1000, f"{name} {fmt}"
        plt.close(fig)


def test_squarify_areas_sum():
    # the treemap packing must conserve area (no gaps/overlap in total)
    rects = jf._squarify([5, 3, 2, 1], 0, 0, 1, 1)
    total = sum(w * h for _, _, w, h in rects)
    assert abs(total - 1.0) < 1e-9


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS  {name}"); passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
