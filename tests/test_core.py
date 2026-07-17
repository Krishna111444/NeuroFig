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


def test_clean_columns_handles_mess():
    import pandas as pd
    df = pd.DataFrame([[1, 2, 3]], columns=[" Group ", "Value\n(mV)", "Unnamed: 2"])
    out = nc.clean_columns(df)
    assert list(out.columns) == ["Group", "Value (mV)", "column_3"]
    # duplicate names get de-duplicated, not silently merged
    df2 = pd.DataFrame([[1, 2]], columns=["x", "x"])
    assert list(nc.clean_columns(df2).columns) == ["x", "x_1"]


def test_wide_to_long_roundtrip():
    import pandas as pd
    wide = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [4.0, 5.0, None]})
    long = nc.wide_to_long(wide, ["A", "B"], group_name="grp", value_name="val")
    assert set(long["grp"]) == {"A", "B"}
    assert len(long) == 5  # the NaN was dropped
    assert "val" in long.columns


def test_nonparametric_forced_two_and_three_groups():
    df = nc.load_table(SAMPLE)
    # 3 groups forced non-parametric -> Kruskal-Wallis + Dunn
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct",
                                  order=["eYFP", "ChR2", "ChR2+Antagonist"],
                                  method="nonparametric")
    assert res.test_name == "Kruskal-Wallis"
    assert res.method == "nonparametric"
    assert res.pairwise and all("p_text" in pr and pr["p_text"] != "0" for pr in res.pairwise)
    # clearly separated groups (eYFP ~17 vs ChR2 ~34) must come out significant
    sig = {(pr["g1"], pr["g2"]): pr["stars"] for pr in res.pairwise}
    key = ("ChR2", "eYFP") if ("ChR2", "eYFP") in sig else ("eYFP", "ChR2")
    assert sig[key] != "ns"

    # 2 groups forced non-parametric -> Mann-Whitney
    df2 = df[df["group"].isin(["eYFP", "ChR2"])]
    res2 = nc.run_group_comparison(df2, "group", "open_arm_time_pct",
                                   order=["eYFP", "ChR2"], method="nonparametric")
    assert res2.test_name == "Mann-Whitney U"


def test_auto_uses_parametric_on_normal_data():
    df = nc.load_table(SAMPLE)
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct",
                                  order=["eYFP", "ChR2", "ChR2+Antagonist"], method="auto")
    # the sample is drawn from normals -> auto should stay parametric
    assert res.method == "parametric"
    assert res.test_name == "One-way ANOVA"


def test_dunn_matches_kruskal_direction():
    # Construct 3 well-separated groups: KW must be significant AND every
    # extreme pair (low vs high) must be significant under Dunn.
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(0)
    rows = []
    for g, mu in [("low", 0), ("mid", 10), ("high", 20)]:
        for v in rng.normal(mu, 1.0, 15):
            rows.append({"g": g, "y": float(v)})
    df = pd.DataFrame(rows)
    res = nc.run_group_comparison(df, "g", "y",
                                  order=["low", "mid", "high"], method="nonparametric")
    assert res.omnibus_p < 0.05
    sig = {frozenset((pr["g1"], pr["g2"])): pr["stars"] for pr in res.pairwise}
    assert sig[frozenset(("low", "high"))] != "ns"


def test_paired_comparison_and_figure():
    import numpy as np
    before = np.array([10.0, 12.0, 9.0, 11.0, 13.0, 10.5, 12.5, 9.5])
    after = before + np.array([2.0, 1.5, 2.5, 1.8, 2.2, 1.9, 2.1, 2.3])  # consistent rise
    res = nc.run_paired_comparison(before, after, method="auto")
    assert res.test_name in ("Paired t-test", "Wilcoxon signed-rank")
    assert res.pairwise[0]["stars"] != "ns"  # consistent within-subject change
    fig = nc.make_paired_figure(before, after, stat=res, ylabel="signal")
    assert len(nc.figure_to_bytes(fig, "pdf")) > 1000


def test_box_and_violin_export():
    df = nc.load_table(SAMPLE)
    order = ["eYFP", "ChR2", "ChR2+Antagonist"]
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct",
                                  order=order, method="nonparametric")
    for kind in ("box", "violin"):
        fig = nc.make_box_figure(df, "group", "open_arm_time_pct", order=order,
                                 kind=kind, stat=res, ylabel="Open-arm time (%)")
        assert len(nc.figure_to_bytes(fig, "svg")) > 1000


def test_timecourse_auc_exact():
    import numpy as np
    # flat trace of height 2 over a width-4 window -> area exactly 8
    t = np.linspace(0, 4, 50)
    flat = np.full((1, 50), 2.0)
    auc = nc.timecourse_auc(t, flat, (0, 4))
    assert abs(float(auc[0]) - 8.0) < 1e-9
    # too-narrow window is rejected, not silently wrong
    try:
        nc.timecourse_auc(t, flat, (0, 0))
        assert False, "expected ValueError for degenerate window"
    except ValueError:
        pass


def test_timecourse_stats_and_figure():
    import numpy as np
    rng = np.random.default_rng(5)
    t = np.linspace(-1, 3, 80)
    A = rng.normal(0, 0.2, (12, t.size)) + (2.0 * (t >= 0))
    B = rng.normal(0, 0.2, (12, t.size)) + (0.5 * (t >= 0))
    mean, sem = nc.timecourse_stats(A, error="sem")
    assert mean.shape == t.shape and sem.shape == t.shape
    assert np.all(sem >= 0)
    # SEM must be smaller than SD for n>1
    _, sd = nc.timecourse_stats(A, error="sd")
    assert np.nanmean(sem) < np.nanmean(sd)
    fig = nc.make_timecourse_figure(t, {"A": A, "B": B}, event_window=(0, 1),
                                    event_label="stim")
    assert len(nc.figure_to_bytes(fig, "pdf")) > 1000


def test_timecourse_from_wide_shape():
    import numpy as np
    import pandas as pd
    df = pd.DataFrame({"t": [0.0, 1.0, 2.0], "s1": [1, 2, 3], "s2": [4, 5, 6]})
    t, m = nc.timecourse_from_wide(df, "t", ["s1", "s2"])
    assert t.tolist() == [0.0, 1.0, 2.0]
    assert m.shape == (2, 3)  # 2 subjects x 3 timepoints


def test_journal_width_is_exact():
    df = nc.load_table(SAMPLE)
    order = ["eYFP", "ChR2", "ChR2+Antagonist"]
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct", order=order)
    fig = nc.make_group_figure(df, "group", "open_arm_time_pct", order=order,
                               stat=res, preset="Nature", width_mm=89.0)
    w_in = fig.get_size_inches()[0]
    assert abs(w_in - 89.0 / 25.4) < 1e-6  # exact to the micrometre


def test_journal_width_lookup_and_fallback():
    assert nc.journal_width_mm("Nature", "single") == 89.0
    assert nc.journal_width_mm("Nature", "double") == 183.0
    assert nc.journal_width_mm("Science", "single") == 55.0
    # unknown journal falls back to Default, does not raise
    assert nc.journal_width_mm("MadeUpJournal") == nc.journal_width_mm("Default")


def test_preset_sets_font_size():
    import matplotlib
    nc.apply_preset("Nature")
    assert matplotlib.rcParams["font.size"] == 7
    nc.apply_preset("Default")
    assert matplotlib.rcParams["font.size"] == 8
    # DejaVu Sans is always the final fallback so text never renders as boxes
    nc.apply_preset("Cell")
    assert "DejaVu Sans" in matplotlib.rcParams["font.sans-serif"]
    nc.apply_preset("Default")  # reset for other tests


def test_license_valid_expired_tampered():
    import licensing as L
    secret, now = "seller-secret", 1_700_000_000
    key = L.make_license("buyer@lab.edu", 30, secret, issued_at=now)
    ok, reason, info = L.verify_license(key, secret, now=now)
    assert ok and reason == "ok" and info["e"] == "buyer@lab.edu"
    # expired
    assert L.verify_license(key, secret, now=now + 31 * 86400)[0] is False
    # tampered signature / wrong secret cannot forge
    assert L.verify_license(key, "not-the-secret", now=now)[0] is False
    body, sig = key.split(".")
    assert L.verify_license(body + ".AAAA", secret, now=now)[0] is False
    # garbage never crashes
    assert L.verify_license("junk", secret)[0] is False
    assert L.verify_license("", secret)[0] is False


def test_watermark_differs_from_clean():
    df = nc.load_table(SAMPLE)
    order = ["eYFP", "ChR2", "ChR2+Antagonist"]
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct", order=order)
    fig = nc.make_group_figure(df, "group", "open_arm_time_pct", order=order, stat=res)
    wm = nc.watermarked_png(fig)
    clean = nc.figure_to_bytes(fig, "png")
    assert len(wm) > 1000 and wm != clean


def test_robustness_guards():
    import numpy as np
    import pandas as pd
    df = nc.load_table(SAMPLE)
    # single group -> clear error, not an opaque scipy crash
    try:
        nc.run_group_comparison(df, "group", "open_arm_time_pct", order=["eYFP"])
        assert False, "expected ValueError for a single group"
    except ValueError:
        pass
    # an empty group is dropped (with a warning), comparison still runs
    res = nc.run_group_comparison(df, "group", "open_arm_time_pct",
                                  order=["eYFP", "ChR2", "GhostGroup"])
    assert any("no data" in w for w in res.warnings)
    assert res.test_name  # a real test still ran on the 2 usable groups
    # paired with <2 pairs -> clear error
    try:
        nc.run_paired_comparison(np.array([1.0]), np.array([2.0]))
        assert False, "expected ValueError for <2 pairs"
    except ValueError:
        pass


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
