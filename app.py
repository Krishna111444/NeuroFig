"""
app.py — the web UI (Streamlit) that wraps neurofig_core.

Run locally:   streamlit run app.py
This is the MVP: upload -> confirm columns -> get figure + stats -> download.
The paywall, accounts, and the plain-English/AI layer are added on top later
(see README.md "Build roadmap").
"""
import streamlit as st
import neurofig_core as nc

st.set_page_config(page_title="NeuroFig — data to publication figure", layout="wide")

st.title("NeuroFig")
st.caption("Upload your data → get a journal-ready figure with the correct statistics. "
           "Colourblind-safe, vector export, no coding.")

# ---- 1. upload -------------------------------------------------------------
up = st.file_uploader("Upload a CSV or Excel file (long format: one column of groups, "
                      "one column of values)", type=["csv", "xlsx", "xls"])

# a built-in sample so first-time visitors can try it with zero friction
use_sample = st.checkbox("No data handy? Use a sample neuroscience dataset", value=up is None)

if up is not None:
    df = nc.load_table(up)
elif use_sample:
    import pandas as pd, numpy as np
    rng = np.random.default_rng(42)
    rows = []
    for g, mu, sd in [("eYFP", 18, 5.5), ("ChR2", 34, 7), ("ChR2+Antagonist", 21, 6)]:
        for _ in range(12):
            rows.append({"group": g, "open_arm_time_pct": round(float(np.clip(rng.normal(mu, sd), 0, 100)), 2)})
    df = pd.DataFrame(rows)
else:
    st.stop()

st.subheader("Your data")
st.dataframe(df.head(20), use_container_width=True)

# ---- 2. confirm columns (pre-filled by the inference engine) --------------
guess = nc.infer_columns(df)
c1, c2 = st.columns(2)
group_col = c1.selectbox("Grouping column", options=list(df.columns),
                         index=list(df.columns).index(guess.group_col) if guess.group_col in df.columns else 0)
value_options = [c for c in df.columns if c != group_col]
default_val = guess.value_col if guess.value_col in value_options else value_options[0]
value_col = c2.selectbox("Value column (numeric)", options=value_options,
                         index=value_options.index(default_val))

levels = list(dict.fromkeys(df[group_col].dropna().tolist()))
order = st.multiselect("Group order (drag to reorder / deselect to drop)", options=levels, default=levels)
ylabel = st.text_input("Y-axis label", value=value_col.replace("_", " "))
title = st.text_input("Panel title (optional)", value="")

o1, o2 = st.columns(2)
plot_type = o1.selectbox("Plot type", ["Grouped bar (mean±SEM)", "Box plot", "Violin plot"])
test_choice = o2.selectbox(
    "Statistical test",
    ["Automatic (recommended)", "Force parametric (t-test / ANOVA)",
     "Force non-parametric (rank-based)"],
    help="Automatic picks parametric vs non-parametric from a normality check, "
         "and always tells you which test it used.")
_method = {"Automatic (recommended)": "auto",
           "Force parametric (t-test / ANOVA)": "parametric",
           "Force non-parametric (rank-based)": "nonparametric"}[test_choice]

# ---- 3. run + render -------------------------------------------------------
if st.button("Generate figure", type="primary") and order:
    stat = nc.run_group_comparison(df, group_col, value_col, order=order, method=_method)
    for w in stat.warnings:
        st.warning(w)
    if plot_type == "Box plot":
        fig = nc.make_box_figure(df, group_col, value_col, order=order, kind="box",
                                 stat=stat, ylabel=ylabel, title=title)
    elif plot_type == "Violin plot":
        fig = nc.make_box_figure(df, group_col, value_col, order=order, kind="violin",
                                 stat=stat, ylabel=ylabel, title=title)
    else:
        fig = nc.make_group_figure(df, group_col, value_col, order=order,
                                   stat=stat, ylabel=ylabel, title=title)
    st.pyplot(fig, use_container_width=False)

    _posthoc = f" · post-hoc: {stat.posthoc_name}" if stat.posthoc_name and len(order) > 2 else ""
    st.success(f"Test used: {stat.test_name} ({stat.method}){_posthoc}  —  {stat.summary}")
    if stat.pairwise:
        st.write("**Pairwise comparisons**")
        st.table([{"Group 1": p["g1"], "Group 2": p["g2"],
                   "p-value": p["p_text"], "Significance": p["stars"]} for p in stat.pairwise])

    # ---- 4. download (vector = the real journal deliverable) --------------
    d1, d2, d3 = st.columns(3)
    d1.download_button("Download PNG", nc.figure_to_bytes(fig, "png"), "figure.png", "image/png")
    d2.download_button("Download PDF (vector)", nc.figure_to_bytes(fig, "pdf"), "figure.pdf", "application/pdf")
    d3.download_button("Download SVG (vector)", nc.figure_to_bytes(fig, "svg"), "figure.svg", "image/svg+xml")

    # ---- 5. where the paywall goes ---------------------------------------
    st.divider()
    st.info("💡 Build note: this is where a paywall gates the high-res / no-watermark "
            "download. Free = watermarked preview; paid = clean vector files.")
