"""
neurofig_core.py
----------------
The deterministic engine of the product. No UI, no AI — just:
  load data -> understand its shape -> run the correct statistics -> draw a
  publication-ready figure.

Keeping this logic separate from the UI (app.py) means it is testable on its own
and easy for Claude Code to extend one function at a time.
"""
from __future__ import annotations
import io
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

try:
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    _HAS_SM = True
except Exception:  # statsmodels optional; ANOVA still works without post-hoc
    _HAS_SM = False

# ---------------------------------------------------------------- palette / style
# Okabe-Ito colourblind-safe palette (the scientific-figure gold standard).
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#7F7F7F"]

JOURNAL_STYLE = {
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 400, "figure.dpi": 130, "pdf.fonttype": 42, "ps.fonttype": 42,
}


def apply_style():
    mpl.rcParams.update(JOURNAL_STYLE)


# ---------------------------------------------------------------- data understanding
@dataclass
class ColumnGuess:
    """What the tool infers about an uploaded table before asking the user."""
    group_col: str | None
    value_col: str | None
    numeric_cols: list[str] = field(default_factory=list)
    categorical_cols: list[str] = field(default_factory=list)


def load_table(file_or_path) -> pd.DataFrame:
    """Accept a path, bytes, or file-like object; read CSV or Excel."""
    if hasattr(file_or_path, "read"):
        name = getattr(file_or_path, "name", "")
        if name.lower().endswith((".xlsx", ".xls")):
            return pd.read_excel(file_or_path)
        return pd.read_csv(file_or_path)
    if str(file_or_path).lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(file_or_path)
    return pd.read_csv(file_or_path)


def infer_columns(df: pd.DataFrame) -> ColumnGuess:
    """Best-effort guess of which column is the grouping factor and which the measure.

    Heuristic: a good grouping column is non-numeric (or low-cardinality) with a
    handful of repeated levels; a good value column is numeric with many values.
    This is deliberately simple — the UI always lets the user override.
    """
    numeric = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    categorical = [c for c in df.columns if c not in numeric]

    # group: prefer a categorical column with 2..8 unique levels
    group_col = None
    best_levels = None
    for c in categorical + numeric:
        nun = df[c].nunique(dropna=True)
        if 2 <= nun <= 8:
            if best_levels is None or nun < best_levels:
                group_col, best_levels = c, nun

    # value: the numeric column with the most distinct values that isn't the group
    value_col = None
    best_card = -1
    for c in numeric:
        if c == group_col:
            continue
        card = df[c].nunique(dropna=True)
        if card > best_card:
            value_col, best_card = c, card

    return ColumnGuess(group_col=group_col, value_col=value_col,
                       numeric_cols=numeric, categorical_cols=categorical)


# ---------------------------------------------------------------- statistics
@dataclass
class StatResult:
    test_name: str
    summary: str                      # human-readable one-liner
    omnibus_p: float | None
    pairwise: list[dict] = field(default_factory=list)  # {g1,g2,p,stars}
    warnings: list[str] = field(default_factory=list)


def _stars(p: float) -> str:
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def fmt_p(p: float) -> str:
    """Journal-standard p-value text. Never prints 'p = 0'.

    Statistical engines (incl. Tukey HSD) clamp extremely small p-values to
    exactly 0.0; reporting that verbatim is scientifically invalid and a
    manuscript red flag. Convention: values below 0.001 are shown as '< 0.001'.
    """
    if p is None or p != p:  # None or NaN
        return "n/a"
    if p < 0.001:
        return "< 0.001"
    if p >= 0.9995:
        return "> 0.999"
    return f"{p:.3f}"


def run_group_comparison(df: pd.DataFrame, group_col: str, value_col: str,
                         order: list[str] | None = None) -> StatResult:
    """Pick and run the appropriate test based on the number of groups.

    2 groups   -> Welch's t-test (does not assume equal variance).
    3+ groups  -> one-way ANOVA, with Tukey HSD post-hoc if statsmodels present.
    Also runs a Shapiro normality check and warns (does not block) on violation.
    """
    d = df[[group_col, value_col]].dropna()
    levels = order or list(pd.unique(d[group_col]))
    samples = [d.loc[d[group_col] == g, value_col].values for g in levels]
    warnings: list[str] = []

    for g, s in zip(levels, samples):
        if len(s) < 3:
            warnings.append(f"Group '{g}' has n={len(s)} (<3) — results unreliable.")
        elif len(s) >= 3:
            try:
                if stats.shapiro(s).pvalue < 0.05:
                    warnings.append(f"Group '{g}' may be non-normal (Shapiro p<0.05); "
                                    f"consider a non-parametric test.")
            except Exception:
                pass

    if len(levels) == 2:
        t, p = stats.ttest_ind(samples[0], samples[1], equal_var=False)
        return StatResult(
            test_name="Welch's t-test",
            summary=f"Welch's t-test: t={t:.2f}, p {fmt_p(p)}",
            omnibus_p=float(p),
            pairwise=[{"g1": levels[0], "g2": levels[1], "p": float(p),
                       "p_text": fmt_p(p), "stars": _stars(p)}],
            warnings=warnings,
        )

    F, p = stats.f_oneway(*samples)
    res = StatResult(
        test_name="One-way ANOVA",
        summary=f"one-way ANOVA: F({len(levels)-1},{len(d)-len(levels)})={F:.2f}, p {fmt_p(p)}",
        omnibus_p=float(p), warnings=warnings,
    )
    if _HAS_SM:
        tuk = pairwise_tukeyhsd(d[value_col], d[group_col])
        tbl = pd.DataFrame(tuk._results_table.data[1:], columns=tuk._results_table.data[0])
        for _, r in tbl.iterrows():
            pv = float(r["p-adj"])
            res.pairwise.append({"g1": str(r["group1"]), "g2": str(r["group2"]),
                                 "p": pv, "p_text": fmt_p(pv), "stars": _stars(pv)})
    return res


# ---------------------------------------------------------------- figure
def make_group_figure(df: pd.DataFrame, group_col: str, value_col: str,
                      order: list[str] | None = None,
                      stat: StatResult | None = None,
                      ylabel: str | None = None, title: str = "") -> plt.Figure:
    """Bar (mean±SEM) + individual points + significance brackets. Returns a Figure."""
    apply_style()
    d = df[[group_col, value_col]].dropna()
    levels = order or list(pd.unique(d[group_col]))
    colors = {g: OKABE_ITO[i % len(OKABE_ITO)] for i, g in enumerate(levels)}

    fig, ax = plt.subplots(figsize=(0.9 * len(levels) + 1.2, 3.0))
    means = []
    for i, g in enumerate(levels):
        v = d.loc[d[group_col] == g, value_col].values
        m = float(np.mean(v)); sem = float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
        means.append(m)
        ax.bar(i, m, width=0.6, color=colors[g], alpha=0.55, edgecolor=colors[g], linewidth=1.2, zorder=1)
        ax.errorbar(i, m, yerr=sem, color="black", lw=1.0, capsize=3, zorder=3)
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.28
        ax.scatter(i + jit, v, s=14, facecolor="white", edgecolor=colors[g], linewidth=1.0, zorder=4)

    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=0)
    ax.set_ylabel(ylabel or value_col.replace("_", " "))
    if title:
        ax.set_title(title, loc="left", fontweight="bold")

    # significance brackets for adjacent, significant pairs
    if stat and stat.pairwise:
        top = max(means) * 1.15 if means else 1.0
        step = (max(means) * 0.12) if means else 0.1
        idx = {g: i for i, g in enumerate(levels)}
        drawn = 0
        for pr in stat.pairwise:
            if pr["stars"] == "ns":
                continue
            if pr["g1"] not in idx or pr["g2"] not in idx:
                continue
            x1, x2 = sorted((idx[pr["g1"]], idx[pr["g2"]]))
            y = top + drawn * step
            ax.plot([x1, x1, x2, x2], [y, y + step * 0.3, y + step * 0.3, y], lw=0.9, color="black")
            ax.text((x1 + x2) / 2, y + step * 0.3, pr["stars"], ha="center", va="bottom", fontsize=9)
            drawn += 1
        ax.set_ylim(0, top + (drawn + 1) * step)

    if stat:
        fig.text(0.5, -0.04, stat.summary + "   (bars mean±SEM; points = observations)",
                 ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig


def figure_to_bytes(fig: plt.Figure, fmt: str = "png") -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.read()
