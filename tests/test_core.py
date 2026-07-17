"""
Portable smoke + correctness tests for neurofig_core.
Run from the project root:  python -m tests.test_core   (or)  pytest -q

No hardcoded paths — uses the bundled sample_data.csv, so it runs on any OS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import neurofig_core as nc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = os.path.join(HERE, "..", "sample_data.csv")


def test_fmt_p_never_zero():
    # The credibility guarantee: no p-value ever renders as "0".
    assert nc.fmt_p(0.0) == "< 0.001"
    assert nc.fmt_p(1e-15) == "< 0.001"
    assert nc.fmt_p(0.0409) == "0.041"
    assert nc.fmt_p(1.0) == "> 0.999"
    assert nc.fmt_p(float("nan")) == "n/a"


def test_infer_and_three_group_anova():
    df = nc.load_table(SAMPLE)
    g = nc.infer_columns(df)
    assert g.group_col == "group"
    assert g.value_col == "open_arm_time_pct"

    order = ["eYFP", "ChR2", "ChR2+Antagonist"]
    res = nc.run_group_comparison(df, g.group_col, g.value_col, order=order)
    assert res.test_name == "One-way ANOVA"
    assert res.pairwise, "expected Tukey pairwise results"
    # every pairwise row must carry a human-safe p_text, never a bare '0'
    for pr in res.pairwise:
        assert pr["p_text"] and pr["p_text"] != "0"


def test_two_group_autoselects_ttest():
    df = nc.load_table(SAMPLE)
    df2 = df[df["group"].isin(["eYFP", "ChR2"])]
    res = nc.run_group_comparison(df2, "group", "open_arm_time_pct",
                                  order=["eYFP", "ChR2"])
    assert res.test_name == "Welch's t-test"


def test_figure_exports_all_formats():
    df = nc.load_table(SAMPLE)
    order = ["eYFP", "ChR2", "ChR2+Antagonist"]
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct", order=order)
    fig = nc.make_group_figure(df, "group", "open_arm_time_pct", order=order,
                               stat=res, ylabel="Open-arm time (%)")
    for fmt in ("png", "pdf", "svg"):
        b = nc.figure_to_bytes(fig, fmt)
        assert len(b) > 1000, f"{fmt} export looks empty"


def _run_all():
    passed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS  {name}")
            passed += 1
    print(f"\n{passed} tests passed.")


if __name__ == "__main__":
    _run_all()
