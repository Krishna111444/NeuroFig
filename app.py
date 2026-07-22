"""
app.py — the NeuroFig web UI (Streamlit).

Tabs:
  • Figures & statistics — group comparison OR ΔF/F time-course (with a
    cluster-based permutation test); upload or type data; edit the figure.
  • PDB → PDBQT — clean a structure and convert it for docking.
  • Calculators — deterministic bench calculations across the life sciences.

Run locally:  streamlit run app.py
"""
import os
import tempfile

import numpy as np
import pandas as pd
import streamlit as st

import cluster_stats as cs
import doseresponse as dr
import journal_figures as jf
import labcalc as lc
import licensing
import neurofig_core as nc
import pdbqt_tools as pq
import standardcurve as sc

st.set_page_config(page_title="NeuroFig — data to publication figure", layout="wide")


def _secret(name: str, default: str = "") -> str:
    try:
        if name in st.secrets:
            return str(st.secrets[name])
    except Exception:
        pass
    return os.environ.get(name, default)


_LICENSE_SECRET = _secret("NEUROFIG_LICENSE_SECRET")
_BUY_URL = _secret("NEUROFIG_BUY_URL")

with st.sidebar:
    st.header("Clean exports")
    st.write("Free preview is watermarked. Unlock clean vector (PDF/SVG) downloads:")
    _key = st.text_input("License key", type="password", placeholder="paste your key")
    unlocked, _reason = False, ""
    if _key:
        if not _LICENSE_SECRET:
            _reason = "Licensing not configured on this server."
        else:
            unlocked, _reason, _ = licensing.verify_license(_key, _LICENSE_SECRET)
        st.success("Unlocked — clean exports enabled.") if unlocked \
            else st.error(f"Key rejected: {_reason}")
    if not unlocked and _BUY_URL:
        st.link_button("Buy a license", _BUY_URL)
    st.caption("Privacy: your data is processed to make the figure and is **not** used "
               "for training and **not** stored after your session.")

st.title("NeuroFig")
st.caption("Data → journal-ready figure with the correct statistics. "
           "Colourblind-safe, vector export, no coding.")


def render_downloads(fig, width_mm, preset, columns):
    """Gated download UI: watermarked free preview vs exact-mm clean vectors."""
    st.divider()
    if unlocked:
        st.write(f"**Clean download** — sized to {width_mm:.0f} mm ({preset}, {columns} column)")
        d1, d2, d3 = st.columns(3)
        d1.download_button("PNG (400 dpi)", nc.figure_to_bytes(fig, "png", exact=True),
                           "figure.png", "image/png")
        d2.download_button("PDF (vector)", nc.figure_to_bytes(fig, "pdf", exact=True),
                           "figure.pdf", "application/pdf")
        d3.download_button("SVG (vector)", nc.figure_to_bytes(fig, "svg", exact=True),
                           "figure.svg", "image/svg+xml")
    else:
        st.write("**Free preview** (watermarked). Unlock clean vector exports in the sidebar.")
        st.download_button("Download watermarked PNG",
                           nc.watermarked_png(fig, "NeuroFig — preview"),
                           "figure_preview.png", "image/png")
        if _BUY_URL:
            st.link_button("Get clean PDF/SVG — buy a license", _BUY_URL)


tab_fig, tab_gallery, tab_pdbqt, tab_calc = st.tabs(
    ["📊 Figures & statistics", "🖼️ Figure gallery", "🧬 PDB → PDBQT", "🧮 Calculators"])


# ============================================================ FIGURES & STATS
with tab_fig:
    analysis = st.radio("Analysis type",
                        ["Group comparison", "Time-course (ΔF/F)",
                         "Dose-response (IC50/EC50)", "Standard curve"], horizontal=True)
    j1, j2 = st.columns(2)
    preset = j1.selectbox("Journal style", list(nc.JOURNAL_PRESETS.keys()), key="preset")
    columns = j2.selectbox("Column width", ["single", "double"], key="cols")
    width_mm = nc.journal_width_mm(preset, columns)

    # ----------------------------------------------------- group comparison
    if analysis == "Group comparison":
        st.subheader("1 · Your data")
        mode = st.radio("Provide data by", ["Upload a file", "Type / paste data"],
                        horizontal=True)
        df = None
        if mode == "Upload a file":
            up = st.file_uploader("CSV/Excel (long: one group column, one value column)",
                                  type=["csv", "xlsx", "xls"])
            use_sample = st.checkbox("Use a sample dataset", value=up is None)
            if up is not None:
                df = nc.load_table(up)
            elif use_sample:
                rng = np.random.default_rng(42)
                rows = []
                for g, mu, sd in [("eYFP", 18, 5.5), ("ChR2", 34, 7), ("ChR2+Antagonist", 21, 6)]:
                    for _ in range(12):
                        rows.append({"group": g,
                                     "open_arm_time_pct": round(float(np.clip(rng.normal(mu, sd), 0, 100)), 2)})
                df = pd.DataFrame(rows)
        else:
            c1, c2 = st.columns(2)
            n_groups = c1.number_input("Number of groups", 2, 12, 3, 1)
            names_raw = c2.text_input("Group names (comma-separated)",
                                      ", ".join(f"Group {i+1}" for i in range(int(n_groups))))
            groups = [g.strip() for g in names_raw.split(",") if g.strip()][:int(n_groups)]
            if len(groups) >= 2:
                edited = st.data_editor(nc.make_entry_template(groups, 10),
                                        num_rows="dynamic", use_container_width=True, key="grid")
                long = nc.entered_wide_to_long(edited, group_cols=groups)
                df = long if not long.empty else None

        if df is None or df.empty:
            st.info("Provide some data above to continue.")
            st.stop()
        st.dataframe(df.head(20), use_container_width=True)

        st.subheader("2 · Comparison & style")
        guess = nc.infer_columns(df)
        a1, a2 = st.columns(2)
        group_col = a1.selectbox("Grouping column", list(df.columns),
                                 index=list(df.columns).index(guess.group_col)
                                 if guess.group_col in df.columns else 0)
        value_options = [c for c in df.columns if c != group_col]
        default_val = guess.value_col if guess.value_col in value_options else value_options[0]
        value_col = a2.selectbox("Value column", value_options,
                                 index=value_options.index(default_val))
        levels = list(dict.fromkeys(df[group_col].dropna().tolist()))
        order = st.multiselect("Group order (deselect to drop)", levels, default=levels)

        o1, o2 = st.columns(2)
        plot_type = o1.selectbox("Plot type", ["Grouped bar (mean±SEM)", "Box plot", "Violin plot"])
        test_choice = o2.selectbox("Statistical test",
                                   ["Automatic (recommended)", "Force parametric (t-test / ANOVA)",
                                    "Force non-parametric (rank-based)"])
        _method = {"Automatic (recommended)": "auto",
                   "Force parametric (t-test / ANOVA)": "parametric",
                   "Force non-parametric (rank-based)": "nonparametric"}[test_choice]

        with st.expander("✏️  Edit figure text"):
            e1, e2 = st.columns(2)
            title = e1.text_input("Title", "", key="gtitle")
            title_loc = e2.selectbox("Title position", ["left", "center", "right"], key="gloc")
            ylabel = e1.text_input("Y-axis label", value=value_col.replace("_", " "), key="gyl")
            relabel = st.text_input("Rename groups (comma-separated, in order)", "", key="grel")

        if len(order) < 2:
            st.warning("Select at least two groups.")
            st.stop()
        try:
            stat = nc.run_group_comparison(df, group_col, value_col, order=order, method=_method)
        except ValueError as e:
            st.error(str(e)); st.stop()
        for w in stat.warnings:
            st.warning(w)

        _kw = dict(order=order, stat=stat, ylabel=ylabel or None, preset=preset, width_mm=width_mm)
        if plot_type == "Box plot":
            fig = nc.make_box_figure(df, group_col, value_col, kind="box", **_kw)
        elif plot_type == "Violin plot":
            fig = nc.make_box_figure(df, group_col, value_col, kind="violin", **_kw)
        else:
            fig = nc.make_group_figure(df, group_col, value_col, **_kw)

        edits = {}
        if title:
            edits.update(title=title, title_loc=title_loc)
        new_labels = [s.strip() for s in relabel.split(",") if s.strip()]
        if new_labels and len(new_labels) == len(order):
            edits["xtick_labels"] = new_labels
        if edits:
            nc.customize_figure(fig, **edits)

        st.subheader("3 · Result")
        st.pyplot(fig, use_container_width=False)
        _ph = f" · post-hoc: {stat.posthoc_name}" if stat.posthoc_name and len(order) > 2 else ""
        st.success(f"Test: {stat.test_name} ({stat.method}){_ph}  —  {stat.summary}")
        if stat.pairwise:
            st.table([{"Group 1": p["g1"], "Group 2": p["g2"],
                       "p-value": p["p_text"], "Significance": p["stars"]} for p in stat.pairwise])
        render_downloads(fig, width_mm, preset, columns)

    # ----------------------------------------------------- time-course ΔF/F
    elif analysis == "Time-course (ΔF/F)":
        st.subheader("1 · Time-course data (wide: one time column + one column per subject)")
        up = st.file_uploader("CSV/Excel", type=["csv", "xlsx", "xls"], key="tcup")
        use_sample = st.checkbox("Use a sample photometry dataset", value=up is None)
        if up is not None:
            tdf = nc.load_table(up)
        elif use_sample:
            rng = np.random.default_rng(4)
            t = np.linspace(-2, 6, 160)
            cols = {"time": t}
            for i in range(12):
                cols[f"CSplus_{i+1}"] = (3.0 * np.exp(-((t - 1.2) ** 2) / (2 * 0.9 ** 2)) * (t >= 0)
                                        + rng.normal(0, 0.3, t.size))
            for i in range(12):
                cols[f"CSminus_{i+1}"] = (1.2 * np.exp(-((t - 1.2) ** 2) / (2 * 0.9 ** 2)) * (t >= 0)
                                         + rng.normal(0, 0.3, t.size))
            tdf = pd.DataFrame(cols)
        else:
            st.info("Provide time-course data to continue."); st.stop()

        st.dataframe(tdf.head(8), use_container_width=True)
        cols_all = list(tdf.columns)
        time_col = st.selectbox("Time column", cols_all,
                                index=0 if "time" not in cols_all else cols_all.index("time"))
        trace_cols = [c for c in cols_all if c != time_col]
        cA = st.multiselect("Condition A columns", trace_cols,
                            default=[c for c in trace_cols if "plus" in c.lower()] or trace_cols[:len(trace_cols)//2])
        cB = st.multiselect("Condition B columns (optional)", trace_cols,
                            default=[c for c in trace_cols if "minus" in c.lower()])
        n1, n2 = st.columns(2)
        labelA = n1.text_input("Condition A name", "CS+")
        labelB = n2.text_input("Condition B name", "CS−")
        e1, e2, e3 = st.columns(3)
        err = e1.selectbox("Error band", ["sem", "sd"])
        show_event = e2.checkbox("Mark stimulus window", value=True)
        ev0 = e3.number_input("Stimulus start", value=0.0) if show_event else None
        ev1 = e3.number_input("Stimulus end", value=2.0) if show_event else None
        ylabel = st.text_input("Y-axis label", "ΔF/F (%)")
        title = st.text_input("Title", "", key="ttitle")

        if not cA:
            st.warning("Select at least the Condition A columns."); st.stop()

        t, A = nc.timecourse_from_wide(tdf, time_col, cA)
        conditions = {labelA: A}
        B = None
        if cB:
            _, B = nc.timecourse_from_wide(tdf, time_col, cB)
            conditions[labelB] = B

        sig_windows = None
        run_cluster = False
        if B is not None:
            st.markdown("**Compare the two conditions across time** — the correct, "
                        "multiple-comparison-controlled test:")
            run_cluster = st.button("Run cluster-based permutation test", type="primary")

        if run_cluster and B is not None:
            with st.spinner("Permuting…"):
                res = cs.cluster_permutation_test(A, B, n_perm=1000, seed=0)
            sig = [c for c in res["clusters"] if c["sig"]]
            sig_windows = [(float(t[c["start"]]), float(t[c["end"]])) for c in sig]
            if sig:
                spans = ", ".join(f"{t[c['start']]:.2f}–{t[c['end']]:.2f}s (p={c['p']:.3f})" for c in sig)
                st.success(f"Significant difference in: {spans}")
            else:
                st.info("No significant cluster — the conditions don't differ after "
                        "controlling for multiple comparisons across time.")

        ev = (ev0, ev1) if show_event else None
        fig = nc.make_timecourse_figure(t, conditions, ylabel=ylabel, error=err,
                                        event_window=ev, event_label="Stimulus" if show_event else None,
                                        title=title, sig_windows=sig_windows,
                                        preset=preset, width_mm=width_mm)
        st.pyplot(fig, use_container_width=False)
        render_downloads(fig, width_mm, preset, columns)

    # ----------------------------------------------------- dose-response
    elif analysis == "Dose-response (IC50/EC50)":
        st.subheader("1 · Dose-response data (a concentration column and a response column)")
        mode = st.radio("Provide data by", ["Upload a file", "Type / paste data"],
                        horizontal=True, key="drmode")
        ddf = None
        if mode == "Upload a file":
            up = st.file_uploader("CSV/Excel", type=["csv", "xlsx", "xls"], key="drup")
            use_sample = st.checkbox("Use a sample drug-screen dataset", value=up is None, key="drsamp")
            if up is not None:
                ddf = nc.load_table(up)
            elif use_sample:
                rng = np.random.default_rng(5)
                conc = np.repeat(np.logspace(-9, -3, 8), 3)
                resp = 100 / (1 + 10 ** ((np.log10(1e-6) - np.log10(conc)) * 1.0)) \
                    + rng.normal(0, 3, conc.size)
                ddf = pd.DataFrame({"concentration_M": conc, "response_pct": np.round(resp, 2)})
        else:
            ddf = st.data_editor(pd.DataFrame({"concentration": [None] * 8, "response": [None] * 8}),
                                 num_rows="dynamic", use_container_width=True, key="drgrid")
            ddf = ddf.apply(pd.to_numeric, errors="coerce").dropna()
            ddf = ddf if not ddf.empty else None

        if ddf is None or ddf.empty:
            st.info("Provide dose-response data to continue."); st.stop()
        st.dataframe(ddf.head(12), use_container_width=True)

        st.subheader("2 · Fit & style")
        cols_all = list(ddf.columns)
        a1, a2 = st.columns(2)
        conc_col = a1.selectbox("Concentration column", cols_all, index=0)
        resp_col = a2.selectbox("Response column", [c for c in cols_all if c != conc_col],
                                index=0)
        kind_label = st.selectbox("Curve type",
                                  ["Auto-detect", "Inhibitor (IC50)", "Agonist (EC50)"])
        kind = {"Auto-detect": "auto", "Inhibitor (IC50)": "ic50",
                "Agonist (EC50)": "ec50"}[kind_label]
        xlabel = a1.text_input("X-axis label", conc_col.replace("_", " "), key="drxl")
        ylabel = a2.text_input("Y-axis label", resp_col.replace("_", " "), key="dryl")
        title = st.text_input("Title", "", key="drtitle")

        try:
            fit = dr.fit_dose_response(ddf[conc_col].to_numpy(), ddf[resp_col].to_numpy(), kind=kind)
        except (ValueError, RuntimeError) as e:
            st.error(str(e)); st.stop()

        st.subheader("3 · Result")
        m1, m2, m3 = st.columns(3)
        ci = fit.ic50_ci
        m1.metric(fit.label, f"{fit.ic50:.3g}",
                  help=f"95% CI {ci[0]:.3g}–{ci[1]:.3g}" if ci[0] == ci[0] else None)
        m2.metric("Hill slope", f"{fit.hill:.2f}")
        m3.metric("R²", f"{fit.r_squared:.3f}")
        if fit.r_squared < 0.9:
            st.warning("Low R² — the 4-parameter logistic may not fit this data well. "
                       "Check for a full sigmoid (clear top and bottom plateaus).")

        fig = dr.make_dose_response_figure(ddf[conc_col].to_numpy(), ddf[resp_col].to_numpy(),
                                           fit, ylabel=ylabel, xlabel=xlabel, title=title,
                                           preset=preset, width_mm=width_mm)
        st.pyplot(fig, use_container_width=False)
        render_downloads(fig, width_mm, preset, columns)

    # ----------------------------------------------------- standard curve
    else:
        st.subheader("1 · Standards (known concentration + measured signal)")
        model_label = st.selectbox("Assay / curve model",
                                   ["Linear (protein assay: BCA/Bradford)",
                                    "qPCR (Ct vs log quantity)", "4PL (ELISA)"])
        model = {"Linear (protein assay: BCA/Bradford)": "linear",
                 "qPCR (Ct vs log quantity)": "qpcr", "4PL (ELISA)": "4pl"}[model_label]

        smode = st.radio("Provide standards by", ["Upload a file", "Type / paste data"],
                         horizontal=True, key="scmode")
        sdf = None
        if smode == "Upload a file":
            up = st.file_uploader("CSV/Excel", type=["csv", "xlsx", "xls"], key="scup")
            use_sample = st.checkbox("Use a sample standard set", value=up is None, key="scsamp")
            if up is not None:
                sdf = nc.load_table(up)
            elif use_sample:
                rng = np.random.default_rng(2)
                if model == "linear":
                    conc = np.repeat([0, 25, 50, 100, 200, 400, 800], 2).astype(float)
                    sig = 0.002 * conc + 0.05 + rng.normal(0, 0.01, conc.size)
                    sdf = pd.DataFrame({"concentration_ug_mL": conc, "A562": np.round(sig, 4)})
                elif model == "qpcr":
                    conc = np.repeat(np.logspace(1, 7, 7), 3)
                    sig = -3.32 * np.log10(conc) + 40 + rng.normal(0, 0.15, conc.size)
                    sdf = pd.DataFrame({"copies": conc, "Ct": np.round(sig, 2)})
                else:
                    conc = np.logspace(-3, 2, 10)
                    sig = 0.05 + 1.95 / (1 + 10 ** ((0 - np.log10(conc)) * 1.2)) + rng.normal(0, 0.03, conc.size)
                    sdf = pd.DataFrame({"analyte_ng_mL": conc, "OD450": np.round(sig, 4)})
        else:
            sdf = st.data_editor(pd.DataFrame({"concentration": [None] * 8, "signal": [None] * 8}),
                                 num_rows="dynamic", use_container_width=True, key="scgrid")
            sdf = sdf.apply(pd.to_numeric, errors="coerce").dropna()
            sdf = sdf if not sdf.empty else None

        if sdf is None or sdf.empty:
            st.info("Provide standards to continue."); st.stop()
        st.dataframe(sdf.head(12), use_container_width=True)

        st.subheader("2 · Fit & unknowns")
        cols_all = list(sdf.columns)
        a1, a2 = st.columns(2)
        conc_col = a1.selectbox("Concentration column", cols_all, index=0, key="scc")
        sig_col = a2.selectbox("Signal column", [c for c in cols_all if c != conc_col], index=0, key="scs")
        unknown_raw = st.text_input("Unknown sample signals (comma-separated, optional)", "",
                                    help="Their concentrations are read off the fitted curve.")
        xlabel = a1.text_input("X label", conc_col.replace("_", " "), key="scxl")
        ylabel = a2.text_input("Y label", sig_col.replace("_", " "), key="scyl")
        title = st.text_input("Title", "", key="sctitle")

        try:
            fit = sc.fit_standard_curve(sdf[conc_col].to_numpy(), sdf[sig_col].to_numpy(), model)
        except (ValueError, RuntimeError) as e:
            st.error(str(e)); st.stop()

        st.subheader("3 · Result")
        m1, m2 = st.columns(2)
        m1.metric("R²", f"{fit.r_squared:.4f}")
        if fit.efficiency is not None:
            m2.metric("Amplification efficiency", f"{fit.efficiency*100:.0f}%",
                      help="90–110% is the accepted qPCR range.")
        if fit.r_squared < 0.98:
            st.warning("R² below 0.98 — check standards for outliers or the wrong model.")

        unknowns = None
        ysig = [float(s) for s in unknown_raw.split(",") if s.strip()]
        if ysig:
            xs = fit.back_calculate(np.array(ysig))
            unknowns = list(zip([f"s{i+1}" for i in range(len(ysig))], ysig))
            st.write("**Back-calculated unknowns**")
            st.table([{"Signal": y, "Concentration": (f"{x:.4g}" if x == x else "out of range")}
                      for y, x in zip(ysig, xs)])

        fig = sc.make_standard_curve_figure(sdf[conc_col].to_numpy(), sdf[sig_col].to_numpy(),
                                            fit, unknowns=unknowns, xlabel=xlabel, ylabel=ylabel,
                                            title=title, preset=preset, width_mm=width_mm)
        st.pyplot(fig, use_container_width=False)
        render_downloads(fig, width_mm, preset, columns)


# ============================================================ FIGURE GALLERY
with tab_gallery:
    st.subheader("20 journal figure types")
    st.caption("Recurring figures from top-journal papers, on the same engine "
               "(colourblind-safe, journal presets, exact-mm vector export). Each is "
               "shown with sample data — pick one, style it, and download.")
    g1, g2, g3 = st.columns([2, 1, 1])
    choice = g1.selectbox("Figure type", list(jf.FIGURES.keys()))
    gpreset = g2.selectbox("Journal style", list(nc.JOURNAL_PRESETS.keys()), key="gal_preset")
    gcols = g3.selectbox("Column width", ["single", "double"], key="gal_cols")
    gwidth = nc.journal_width_mm(gpreset, gcols)

    if choice in jf.REAL_DATA_FIGURES:
        kind = jf.REAL_DATA_FIGURES[choice]
        sample = {"volcano": jf.sample_volcano_df, "raincloud": jf.sample_raincloud_df,
                  "correlation": jf.sample_correlation_df}[kind]()
        st.info("This figure accepts your own data. Download the template to see the "
                "expected format, or upload a matching file.")
        st.download_button("⬇ Sample template (CSV)", sample.to_csv(index=False),
                           f"{kind}_template.csv", "text/csv")
        src = st.radio("Data", ["Use sample", "Upload your data"], horizontal=True, key="gal_src")
        df = sample
        if src == "Upload your data":
            up = st.file_uploader("CSV/Excel", type=["csv", "xlsx", "xls"], key="gal_up")
            df = nc.load_table(up) if up is not None else None
        if df is None:
            st.stop()
        st.dataframe(df.head(10), use_container_width=True)
        cols = list(df.columns)
        gtitle = st.text_input("Title (optional)", "", key="gal_title")

        fig, warns = None, []
        try:
            if kind == "volcano":
                c1, c2, c3 = st.columns(3)
                fc_col = c1.selectbox("Fold-change column", cols)
                p_col = c2.selectbox("p-value column", [c for c in cols if c != fc_col])
                label_col = c3.selectbox("Label column (optional)", ["(none)"] + cols)
                o1, o2, o3, o4 = st.columns(4)
                fc_is_log2 = o1.selectbox("Fold change is", ["already log2", "linear"]) == "already log2"
                p_is_nlp = o2.selectbox("p-value is", ["raw p", "already -log10"]) == "already -log10"
                fc_thresh = o3.number_input("FC threshold", value=1.0, min_value=0.0)
                p_thresh = o4.number_input("p threshold", value=0.05, min_value=0.0, format="%.4f")
                fig, warns = jf.volcano_from_data(
                    df, fc_col, p_col, fc_is_log2=fc_is_log2, p_is_neglog10=p_is_nlp,
                    fc_thresh=fc_thresh, p_thresh=p_thresh,
                    label_col=None if label_col == "(none)" else label_col,
                    preset=gpreset, width_mm=gwidth)
            elif kind == "raincloud":
                c1, c2 = st.columns(2)
                group_col = c1.selectbox("Group column", cols)
                value_col = c2.selectbox("Value column", [c for c in cols if c != group_col])
                fig, warns = jf.raincloud_from_data(df, group_col, value_col,
                                                    preset=gpreset, width_mm=gwidth)
            else:  # correlation
                value_cols = st.multiselect("Numeric columns to correlate", cols, default=cols)
                method = st.selectbox("Method", ["pearson", "spearman", "kendall"])
                if len(value_cols) < 2:
                    st.warning("Select at least two columns."); st.stop()
                fig, warns = jf.correlation_heatmap_from_data(df, value_cols, method=method,
                                                              preset=gpreset, width_mm=gwidth)
        except (ValueError, RuntimeError) as e:
            st.error(str(e)); st.stop()
        for m in warns:
            st.warning(m)
        if gtitle:
            try:
                nc.customize_figure(fig, title=gtitle)
            except Exception:
                pass
        st.pyplot(fig, use_container_width=False)
        render_downloads(fig, gwidth, gpreset, gcols)
    else:
        gtitle = st.text_input("Title (optional)", "", key="gal_title")
        fig = jf.FIGURES[choice](preset=gpreset, width_mm=gwidth)
        if gtitle:
            try:
                nc.customize_figure(fig, title=gtitle)
            except Exception:
                pass
        st.caption("Shown with sample data. Volcano, raincloud, and correlation heatmap "
                   "also accept your own uploads.")
        st.pyplot(fig, use_container_width=False)
        render_downloads(fig, gwidth, gpreset, gcols)


# ============================================================ PDB → PDBQT
with tab_pdbqt:
    st.subheader("Prepare a structure for AutoDock / Vina")
    st.caption("Clean the PDB (remove waters, pick chains, drop alternate atoms), "
               "then convert to PDBQT with correct charges and atom types.")
    pdb_up = st.file_uploader("Upload a .pdb file", type=["pdb", "ent"])
    if pdb_up is None:
        st.info("Upload a PDB file to begin.")
    else:
        text = pdb_up.getvalue().decode("utf-8", errors="replace")
        s = pq.summarize_pdb(text)
        st.write(f"**Contents:** {s.n_atoms} protein atoms · {s.n_hetatm} heteroatoms "
                 f"({s.n_waters} water) · chains {s.chains or '—'} · "
                 f"hetero {s.hetero_residues or '—'} · {s.n_models} model(s)")
        c1, c2 = st.columns(2)
        keep_chains = c1.multiselect("Keep chains", s.chains, default=s.chains)
        remove_waters = c2.checkbox("Remove waters", value=True)
        remove_hetero = c1.checkbox("Remove other heteroatoms", value=True)
        keep_hetero = c2.multiselect("…except keep these residues", s.hetero_residues)
        drop_altloc = st.checkbox("Drop alternate-location duplicates", value=True)
        cleaned, report = pq.clean_pdb(text, keep_chains=keep_chains or None,
                                       remove_waters=remove_waters, remove_hetero=remove_hetero,
                                       keep_hetero=keep_hetero or None, drop_altloc=drop_altloc)
        st.caption("Removed — " + ", ".join(f"{k}: {v}" for k, v in report.items()))
        st.download_button("Download cleaned PDB", cleaned, "cleaned.pdb", "chemical/x-pdb")

        st.divider()
        conv_mode = st.radio("Convert as", ["receptor", "ligand"], horizontal=True)
        ph = st.number_input("Add hydrogens for pH", 0.0, 14.0, 7.4, 0.1)
        if pq.obabel_available():
            if st.button("Convert to PDBQT", type="primary"):
                try:
                    with tempfile.TemporaryDirectory() as d:
                        ip, op = os.path.join(d, "in.pdb"), os.path.join(d, "out.pdbqt")
                        with open(ip, "w") as fh:
                            fh.write(cleaned)
                        pq.pdb_to_pdbqt(ip, op, mode=conv_mode, ph=ph)
                        with open(op) as fh:
                            result = fh.read()
                    st.success("Converted.")
                    st.download_button("Download PDBQT", result, "structure.pdbqt", "text/plain")
                except Exception as e:
                    st.error(f"Conversion failed: {e}")
        else:
            st.warning("Open Babel isn't installed here, so PDBQT conversion is unavailable "
                       "(cleaning still works). To enable it:")
            st.code("conda create -n neurofig_chem -c conda-forge python=3.11 openbabel -y", "bash")


# ============================================================ CALCULATORS
with tab_calc:
    st.subheader("Bench calculators")
    st.caption("Deterministic, reference-validated. No AI, no guessing — just the right number.")
    cat = st.selectbox("Calculator", [
        "Solution: mass for molarity", "Solution: dilution (C1V1=C2V2)",
        "DNA: concentration from A260", "DNA: ng ↔ pmol", "DNA: copy number (qPCR)",
        "DNA: GC% / reverse-complement / translate",
        "Protein: MW & extinction coefficient", "Spectrophotometry: Beer–Lambert",
        "Buffer: Henderson–Hasselbalch", "Microbiology: CFU/mL",
        "Microbiology: doubling time", "Pharmacology: human-equivalent dose"])

    if cat == "Solution: mass for molarity":
        c1, c2, c3 = st.columns(3)
        conc = c1.number_input("Concentration (mol/L)", value=0.1, format="%.4f")
        vol = c2.number_input("Volume (L)", value=1.0, format="%.4f")
        mw = c3.number_input("Molecular weight (g/mol)", value=58.44, format="%.4f")
        st.metric("Mass needed", f"{lc.mass_for_molarity(conc, vol, mw):.4g} g")

    elif cat == "Solution: dilution (C1V1=C2V2)":
        c1, c2, c3 = st.columns(3)
        sc = c1.number_input("Stock concentration", value=10.0)
        tc = c2.number_input("Target concentration", value=1.0)
        tv = c3.number_input("Target volume", value=10.0)
        try:
            v1 = lc.stock_volume_for_target(sc, tc, tv)
            st.metric("Stock volume", f"{v1:.4g}")
            st.caption(f"Add {tv - v1:.4g} of diluent.")
        except ValueError as e:
            st.error(str(e))

    elif cat == "DNA: concentration from A260":
        c1, c2, c3 = st.columns(3)
        a260 = c1.number_input("A260", value=1.0)
        dil = c2.number_input("Dilution factor", value=1.0)
        a280 = c3.number_input("A280 (optional, for purity)", value=0.0)
        st.metric("dsDNA concentration", f"{lc.dsdna_conc_ng_uL(a260, dil):.4g} ng/µL")
        if a280 > 0:
            st.caption(f"A260/A280 purity ratio: {lc.nucleic_purity(a260, a280):.2f} (pure DNA ≈ 1.8)")

    elif cat == "DNA: ng ↔ pmol":
        c1, c2 = st.columns(2)
        ng = c1.number_input("Mass (ng)", value=100.0)
        length = c2.number_input("Length (bp)", value=1000, step=1)
        st.metric("Amount", f"{lc.dna_ng_to_pmol(ng, int(length)):.4g} pmol")

    elif cat == "DNA: copy number (qPCR)":
        c1, c2 = st.columns(2)
        ng = c1.number_input("Mass (ng)", value=1.0, format="%.4f")
        length = c2.number_input("Length (bp)", value=1000, step=1)
        st.metric("Molecule count", f"{lc.dna_copies(ng, int(length)):.4g} copies")

    elif cat == "DNA: GC% / reverse-complement / translate":
        seq = st.text_area("DNA sequence", "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG").strip()
        if seq:
            try:
                st.write(f"**GC content:** {lc.gc_content(seq)*100:.1f}%")
                st.code("reverse complement:\n" + lc.reverse_complement(seq))
                st.code("translation (frame 1):\n" + lc.translate(seq))
            except ValueError as e:
                st.error(str(e))

    elif cat == "Protein: MW & extinction coefficient":
        seq = st.text_area("Protein sequence (1-letter)", "MKWVTFISLLLLFSSAYS").strip()
        cystines = st.number_input("Disulfide bonds (cystines)", value=0, step=1)
        if seq:
            try:
                st.metric("Average MW", f"{lc.protein_average_mw(seq):,.1f} Da")
                st.metric("ε₂₈₀", f"{lc.extinction_coeff_280(seq, int(cystines)):,} M⁻¹cm⁻¹")
            except ValueError as e:
                st.error(str(e))

    elif cat == "Spectrophotometry: Beer–Lambert":
        c1, c2, c3 = st.columns(3)
        a = c1.number_input("Absorbance", value=0.5)
        eps = c2.number_input("ε (M⁻¹cm⁻¹)", value=6990.0)
        path = c3.number_input("Path length (cm)", value=1.0)
        try:
            st.metric("Concentration", f"{lc.beer_lambert_concentration(a, eps, path):.4g} M")
        except ValueError as e:
            st.error(str(e))

    elif cat == "Buffer: Henderson–Hasselbalch":
        c1, c2, c3 = st.columns(3)
        pka = c1.number_input("pKa", value=4.76)
        base = c2.number_input("[A⁻] (base)", value=1.0)
        acid = c3.number_input("[HA] (acid)", value=1.0)
        try:
            st.metric("pH", f"{lc.henderson_hasselbalch_ph(pka, base, acid):.2f}")
        except ValueError as e:
            st.error(str(e))

    elif cat == "Microbiology: CFU/mL":
        c1, c2, c3 = st.columns(3)
        col = c1.number_input("Colonies counted", value=50, step=1)
        dil = c2.number_input("Dilution factor", value=1e-6, format="%.1e")
        pl = c3.number_input("Volume plated (mL)", value=0.1)
        try:
            st.metric("CFU/mL", f"{lc.cfu_per_ml(int(col), dil, pl):.4g}")
        except ValueError as e:
            st.error(str(e))

    elif cat == "Microbiology: doubling time":
        c1, c2, c3 = st.columns(3)
        n0 = c1.number_input("N₀", value=1.0)
        nt = c2.number_input("Nₜ", value=8.0)
        hrs = c3.number_input("Elapsed time", value=6.0)
        try:
            st.metric("Doubling time", f"{lc.doubling_time(n0, nt, hrs):.4g}")
        except ValueError as e:
            st.error(str(e))

    elif cat == "Pharmacology: human-equivalent dose":
        c1, c2 = st.columns(2)
        dose = c1.number_input("Animal dose (mg/kg)", value=10.0)
        sp = c2.selectbox("Species", [k for k in lc.KM_FACTORS if k != "human"])
        st.metric("Human-equivalent dose", f"{lc.human_equivalent_dose(dose, sp):.4g} mg/kg")
        st.caption("FDA body-surface-area scaling (2005 guidance).")
