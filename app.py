"""
app.py — the NeuroFig web UI (Streamlit) wrapping neurofig_core + pdbqt_tools.

Two tabs:
  • Figures & statistics — upload OR type data → pick a graph/test → the engine
    computes the right statistic and draws a journal-ready figure you can edit.
  • PDB → PDBQT — clean a structure and convert it for AutoDock/Vina docking.

Run locally:  streamlit run app.py
"""
import os
import tempfile

import pandas as pd
import streamlit as st

import licensing
import neurofig_core as nc
import pdbqt_tools as pq

st.set_page_config(page_title="NeuroFig — data to publication figure", layout="wide")


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

tab_fig, tab_pdbqt = st.tabs(["📊 Figures & statistics", "🧬 PDB → PDBQT"])


# ============================================================ FIGURES & STATS
with tab_fig:
    # ---- 1. get data: upload OR type it in -------------------------------
    st.subheader("1 · Your data")
    mode = st.radio("How do you want to provide data?",
                    ["Upload a file", "Type / paste data"], horizontal=True)

    df = None
    if mode == "Upload a file":
        up = st.file_uploader("CSV or Excel (long format: one group column, one value column)",
                              type=["csv", "xlsx", "xls"])
        use_sample = st.checkbox("No data handy? Use a sample neuroscience dataset",
                                 value=up is None)
        if up is not None:
            df = nc.load_table(up)
        elif use_sample:
            import numpy as np
            rng = np.random.default_rng(42)
            rows = []
            for g, mu, sd in [("eYFP", 18, 5.5), ("ChR2", 34, 7), ("ChR2+Antagonist", 21, 6)]:
                for _ in range(12):
                    rows.append({"group": g,
                                 "open_arm_time_pct": round(float(np.clip(rng.normal(mu, sd), 0, 100)), 2)})
            df = pd.DataFrame(rows)
    else:
        st.write("Enter one column per group, then type or paste your numbers:")
        c1, c2 = st.columns(2)
        n_groups = c1.number_input("Number of groups", 2, 12, 3, 1)
        names_raw = c2.text_input("Group names (comma-separated)",
                                  ", ".join(["Group %d" % (i + 1) for i in range(int(n_groups))]))
        groups = [g.strip() for g in names_raw.split(",") if g.strip()][:int(n_groups)]
        if len(groups) >= 2:
            template = nc.make_entry_template(groups, n_rows=10)
            edited = st.data_editor(template, num_rows="dynamic", use_container_width=True,
                                    key="grid")
            long = nc.entered_wide_to_long(edited, group_cols=groups)
            df = long if not long.empty else None
            if df is not None:
                st.caption(f"{len(df)} values entered across {df['group'].nunique()} groups.")

    if df is None or df.empty:
        st.info("Provide some data above to continue.")
        st.stop()

    st.dataframe(df.head(20), use_container_width=True)

    # ---- 2. columns + options -------------------------------------------
    st.subheader("2 · Comparison & style")
    guess = nc.infer_columns(df)
    a1, a2 = st.columns(2)
    group_col = a1.selectbox("Grouping column", list(df.columns),
                             index=list(df.columns).index(guess.group_col)
                             if guess.group_col in df.columns else 0)
    value_options = [c for c in df.columns if c != group_col]
    default_val = guess.value_col if guess.value_col in value_options else value_options[0]
    value_col = a2.selectbox("Value column (numeric)", value_options,
                             index=value_options.index(default_val))

    levels = list(dict.fromkeys(df[group_col].dropna().tolist()))
    order = st.multiselect("Group order (deselect to drop)", levels, default=levels)

    o1, o2 = st.columns(2)
    plot_type = o1.selectbox("Plot type", ["Grouped bar (mean±SEM)", "Box plot", "Violin plot"])
    test_choice = o2.selectbox(
        "Statistical test",
        ["Automatic (recommended)", "Force parametric (t-test / ANOVA)",
         "Force non-parametric (rank-based)"],
        help="Automatic picks parametric vs non-parametric from a normality check "
             "and always tells you which test it used.")
    _method = {"Automatic (recommended)": "auto",
               "Force parametric (t-test / ANOVA)": "parametric",
               "Force non-parametric (rank-based)": "nonparametric"}[test_choice]

    j1, j2 = st.columns(2)
    preset = j1.selectbox("Journal style", list(nc.JOURNAL_PRESETS.keys()))
    columns = j2.selectbox("Column width", ["single", "double"])
    width_mm = nc.journal_width_mm(preset, columns)

    # ---- 3. edit the figure's text --------------------------------------
    with st.expander("✏️  Edit figure text (title, labels, group names)"):
        e1, e2 = st.columns(2)
        title = e1.text_input("Title", "")
        title_loc = e2.selectbox("Title position", ["left", "center", "right"])
        ylabel = e1.text_input("Y-axis label", value=value_col.replace("_", " "))
        xlabel = e2.text_input("X-axis label", "")
        relabel = st.text_input("Rename groups (comma-separated, in the order above)", "",
                                help="These labels are the 'legend' of a bar/box plot.")

    if len(order) < 2:
        st.warning("Select at least two groups to compare.")
        st.stop()

    # ---- 4. compute + render --------------------------------------------
    try:
        stat = nc.run_group_comparison(df, group_col, value_col, order=order, method=_method)
    except ValueError as e:
        st.error(str(e))
        st.stop()
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
    if xlabel:
        edits["xlabel"] = xlabel
    new_labels = [s.strip() for s in relabel.split(",") if s.strip()]
    if new_labels and len(new_labels) == len(order):
        edits["xtick_labels"] = new_labels
    if edits:
        nc.customize_figure(fig, **edits)

    st.subheader("3 · Result")
    st.pyplot(fig, use_container_width=False)

    _posthoc = f" · post-hoc: {stat.posthoc_name}" if stat.posthoc_name and len(order) > 2 else ""
    st.success(f"Test used: {stat.test_name} ({stat.method}){_posthoc}  —  {stat.summary}")
    if stat.pairwise:
        st.write("**Pairwise comparisons**")
        st.table([{"Group 1": p["g1"], "Group 2": p["g2"],
                   "p-value": p["p_text"], "Significance": p["stars"]} for p in stat.pairwise])

    # ---- 5. download — free preview watermarked; clean exports gated -----
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
                 f"hetero residues {s.hetero_residues or '—'} · {s.n_models} model(s)")

        c1, c2 = st.columns(2)
        keep_chains = c1.multiselect("Keep chains", s.chains, default=s.chains)
        remove_waters = c2.checkbox("Remove waters", value=True)
        remove_hetero = c1.checkbox("Remove other heteroatoms (ligands, ions)", value=True)
        keep_hetero = c2.multiselect("…except keep these residues", s.hetero_residues)
        drop_altloc = st.checkbox("Drop alternate-location duplicate atoms", value=True)

        cleaned, report = pq.clean_pdb(
            text, keep_chains=keep_chains or None, remove_waters=remove_waters,
            remove_hetero=remove_hetero, keep_hetero=keep_hetero or None,
            drop_altloc=drop_altloc)
        st.caption("Removed — " + ", ".join(f"{k}: {v}" for k, v in report.items()))
        st.download_button("Download cleaned PDB", cleaned, "cleaned.pdb", "chemical/x-pdb")

        st.divider()
        mode = st.radio("Convert as", ["receptor", "ligand"], horizontal=True)
        ph = st.number_input("Add hydrogens for pH", 0.0, 14.0, 7.4, 0.1)
        if pq.obabel_available():
            if st.button("Convert to PDBQT", type="primary"):
                try:
                    with tempfile.TemporaryDirectory() as d:
                        ip = os.path.join(d, "in.pdb"); op = os.path.join(d, "out.pdbqt")
                        with open(ip, "w") as fh:
                            fh.write(cleaned)
                        pq.pdb_to_pdbqt(ip, op, mode=mode, ph=ph)
                        with open(op) as fh:
                            result = fh.read()
                    st.success("Converted.")
                    st.download_button("Download PDBQT", result, "structure.pdbqt", "text/plain")
                except Exception as e:
                    st.error(f"Conversion failed: {e}")
        else:
            st.warning("Open Babel isn't installed on this server, so PDBQT conversion is "
                       "unavailable here. Cleaning above still works. To enable conversion:")
            st.code("conda create -n neurofig_chem -c conda-forge python=3.11 openbabel -y",
                    language="bash")
