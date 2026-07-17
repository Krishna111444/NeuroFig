"""
app.py — the web UI (Streamlit) that wraps neurofig_core.

Run locally:   streamlit run app.py
This is the MVP: upload -> confirm columns -> get figure + stats -> download.
The paywall, accounts, and the plain-English/AI layer are added on top later
(see README.md "Build roadmap").
"""
import os

import streamlit as st

import licensing
import neurofig_core as nc

st.set_page_config(page_title="NeuroFig — data to publication figure", layout="wide")

st.title("NeuroFig")
st.caption("Upload your data → get a journal-ready figure with the correct statistics. "
           "Colourblind-safe, vector export, no coding.")


def _secret(name: str, default: str = "") -> str:
    """Read config from Streamlit secrets first, then environment."""
    try:
        if name in st.secrets:
            return str(st.secrets[name])
    except Exception:
        pass
    return os.environ.get(name, default)


# ---- paywall: unlock clean exports with a license key ----------------------
_LICENSE_SECRET = _secret("NEUROFIG_LICENSE_SECRET")
_BUY_URL = _secret("NEUROFIG_BUY_URL")  # your Lemon Squeezy / Razorpay payment link

with st.sidebar:
    st.header("Clean exports")
    st.write("Free preview is watermarked. Unlock clean vector (PDF/SVG) downloads:")
    _key = st.text_input("License key", type="password", placeholder="paste your key")
    unlocked, _reason = False, ""
    if _key:
        if not _LICENSE_SECRET:
            _reason = "Licensing not configured on this server."
        else:
            unlocked, _reason, _info = licensing.verify_license(_key, _LICENSE_SECRET)
        st.success("Unlocked — clean exports enabled.") if unlocked else st.error(f"Key rejected: {_reason}")
    if not unlocked and _BUY_URL:
        st.link_button("Buy a license", _BUY_URL)
    st.caption("Privacy: your data is processed to make the figure and is **not** used "
               "for training and **not** stored after your session.")

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

j1, j2 = st.columns(2)
preset = j1.selectbox("Journal style", list(nc.JOURNAL_PRESETS.keys()),
                      help="Sets fonts and the exact figure width for that journal's "
                           "single/double column. Vector downloads are sized to the mm.")
columns = j2.selectbox("Column width", ["single", "double"])
width_mm = nc.journal_width_mm(preset, columns)
st.caption(f"Figure will export at **{width_mm:.0f} mm** wide ({preset}, {columns} column).")

# ---- 3. run + render -------------------------------------------------------
if st.button("Generate figure", type="primary") and order:
    stat = nc.run_group_comparison(df, group_col, value_col, order=order, method=_method)
    for w in stat.warnings:
        st.warning(w)
    _kw = dict(order=order, stat=stat, ylabel=ylabel, title=title,
               preset=preset, width_mm=width_mm)
    if plot_type == "Box plot":
        fig = nc.make_box_figure(df, group_col, value_col, kind="box", **_kw)
    elif plot_type == "Violin plot":
        fig = nc.make_box_figure(df, group_col, value_col, kind="violin", **_kw)
    else:
        fig = nc.make_group_figure(df, group_col, value_col, **_kw)
    st.pyplot(fig, use_container_width=False)

    _posthoc = f" · post-hoc: {stat.posthoc_name}" if stat.posthoc_name and len(order) > 2 else ""
    st.success(f"Test used: {stat.test_name} ({stat.method}){_posthoc}  —  {stat.summary}")
    if stat.pairwise:
        st.write("**Pairwise comparisons**")
        st.table([{"Group 1": p["g1"], "Group 2": p["g2"],
                   "p-value": p["p_text"], "Significance": p["stars"]} for p in stat.pairwise])

    # ---- 4. download — free preview watermarked; clean exports gated -------
    st.divider()
    if unlocked:
        # exact=True preserves the journal column width to the millimetre.
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
