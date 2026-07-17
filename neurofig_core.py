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
    "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "savefig.dpi": 400, "figure.dpi": 130, "pdf.fonttype": 42, "ps.fonttype": 42,
}

# Journal figure specs. Column widths in mm are from each journal's published
# author guidelines; base font size follows their typical label size. Arial/
# Helvetica are the requested faces — matplotlib falls back down the list to the
# always-present DejaVu Sans if they are not installed (no crash, no ugly boxes).
_SANS = ["Arial", "Helvetica", "DejaVu Sans"]
JOURNAL_PRESETS: dict[str, dict] = {
    "Default": {"single_mm": 89,  "double_mm": 183, "font_pt": 8, "family": _SANS},
    "Nature":  {"single_mm": 89,  "double_mm": 183, "font_pt": 7, "family": _SANS},
    "Cell":    {"single_mm": 85,  "double_mm": 174, "font_pt": 7, "family": _SANS},
    "eLife":   {"single_mm": 85,  "double_mm": 170, "font_pt": 8, "family": _SANS},
    "Science": {"single_mm": 55,  "double_mm": 120, "font_pt": 7, "family": _SANS},
}

_ACTIVE_PRESET = "Default"


def _rcparams_for(name: str) -> dict:
    p = JOURNAL_PRESETS.get(name, JOURNAL_PRESETS["Default"])
    rc = dict(JOURNAL_STYLE)
    rc.update({"font.family": "sans-serif", "font.sans-serif": p["family"],
               "font.size": p["font_pt"]})
    return rc


def apply_preset(name: str = "Default") -> None:
    """Select a journal style (fonts + sizes). Width is applied per-figure."""
    global _ACTIVE_PRESET
    _ACTIVE_PRESET = name if name in JOURNAL_PRESETS else "Default"
    mpl.rcParams.update(_rcparams_for(_ACTIVE_PRESET))


def apply_style() -> None:
    """Apply the currently-active preset's rcParams (figure functions call this)."""
    mpl.rcParams.update(_rcparams_for(_ACTIVE_PRESET))


def journal_width_mm(name: str, columns: str = "single") -> float:
    p = JOURNAL_PRESETS.get(name, JOURNAL_PRESETS["Default"])
    return float(p["double_mm"] if columns == "double" else p["single_mm"])


def set_width_mm(fig: plt.Figure, width_mm: float) -> plt.Figure:
    """Resize a figure to an exact width in mm, preserving its aspect ratio."""
    w_in = width_mm / 25.4
    cw, ch = fig.get_size_inches()
    fig.set_size_inches(w_in, w_in * ch / cw)
    return fig


# ---------------------------------------------------------------- data understanding
@dataclass
class ColumnGuess:
    """What the tool infers about an uploaded table before asking the user."""
    group_col: str | None
    value_col: str | None
    numeric_cols: list[str] = field(default_factory=list)
    categorical_cols: list[str] = field(default_factory=list)


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Tidy messy real-world headers: strip whitespace, collapse repeats, dedupe.

    Real uploads have columns like ' Group ', 'Value\n(mV)', or two blank names.
    We normalise names but never touch the data, so downstream logic is stable.
    """
    new_cols, seen = [], {}
    for i, c in enumerate(df.columns):
        name = " ".join(str(c).split()).strip()  # collapse internal/edge whitespace
        if not name or name.lower().startswith("unnamed"):
            name = f"column_{i+1}"
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        new_cols.append(name)
    out = df.copy()
    out.columns = new_cols
    return out


def load_table(file_or_path) -> pd.DataFrame:
    """Accept a path, bytes, or file-like object; read CSV or Excel; tidy headers."""
    if hasattr(file_or_path, "read"):
        name = getattr(file_or_path, "name", "")
        if name.lower().endswith((".xlsx", ".xls")):
            df = pd.read_excel(file_or_path)
        else:
            df = pd.read_csv(file_or_path)
    elif str(file_or_path).lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(file_or_path)
    else:
        df = pd.read_csv(file_or_path)
    return clean_columns(df)


def wide_to_long(df: pd.DataFrame, value_cols: list[str],
                 group_name: str = "group", value_name: str = "value") -> pd.DataFrame:
    """Reshape a wide table (one column per group) into long format.

    Prism users routinely have data as 'one column per condition'. The engine
    works in long format, so this bridges the common wide layout to it.
    """
    long = df[value_cols].melt(var_name=group_name, value_name=value_name)
    return long.dropna(subset=[value_name]).reset_index(drop=True)


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
    pairwise: list[dict] = field(default_factory=list)  # {g1,g2,p,p_text,stars}
    warnings: list[str] = field(default_factory=list)
    method: str = "parametric"        # which family was actually used
    posthoc_name: str = ""            # name of the pairwise/post-hoc test


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


def _check_normality(levels, samples) -> tuple[bool, list[str]]:
    """Per-group Shapiro-Wilk. Returns (all_groups_normal, warnings).

    Shapiro is a flag, not a verdict: it is under-powered at small n and
    over-sensitive at large n. We use it only to decide the *default* family.
    """
    warnings: list[str] = []
    all_normal = True
    for g, s in zip(levels, samples):
        if len(s) < 3:
            warnings.append(f"Group '{g}' has n={len(s)} (<3) — results unreliable.")
            all_normal = False
            continue
        try:
            if stats.shapiro(s).pvalue < 0.05:
                all_normal = False
        except Exception:
            pass
    return all_normal, warnings


def _dunn_test(values: np.ndarray, groups: np.ndarray, levels: list) -> list[dict]:
    """Dunn's test (rank-based post-hoc for Kruskal-Wallis), Holm-adjusted.

    Uses global average ranks with a tie correction, the standard formulation
    (Dunn 1964). Holm-Bonferroni controls family-wise error less conservatively
    than plain Bonferroni. Returns one dict per group pair.
    """
    import itertools
    N = len(values)
    ranks = stats.rankdata(values)  # average ranks handle ties
    _, counts = np.unique(values, return_counts=True)
    tie_sum = float(np.sum(counts.astype(float) ** 3 - counts))
    sigma2 = (N * (N + 1) / 12.0) - (tie_sum / (12.0 * (N - 1))) if N > 1 else 0.0

    n = {lv: int((groups == lv).sum()) for lv in levels}
    mean_rank = {lv: float(ranks[groups == lv].mean()) if n[lv] else float("nan")
                 for lv in levels}

    pairs = list(itertools.combinations(levels, 2))
    raw = []
    for a, b in pairs:
        se = np.sqrt(sigma2 * (1.0 / n[a] + 1.0 / n[b])) if (n[a] and n[b]) else 0.0
        z = (mean_rank[a] - mean_rank[b]) / se if se > 0 else 0.0
        p = 2.0 * stats.norm.sf(abs(z))
        raw.append([a, b, z, min(p, 1.0)])

    # Holm-Bonferroni step-down adjustment, with monotonicity enforced
    m = len(raw)
    order_idx = sorted(range(m), key=lambda i: raw[i][3])
    adj = [1.0] * m
    prev = 0.0
    for rank_i, idx in enumerate(order_idx):
        p_adj = min(1.0, (m - rank_i) * raw[idx][3])
        p_adj = max(p_adj, prev)  # ensure non-decreasing
        prev = p_adj
        adj[idx] = p_adj

    out = []
    for (a, b, z, _), pa in zip(raw, adj):
        out.append({"g1": str(a), "g2": str(b), "p": float(pa),
                    "p_text": fmt_p(pa), "stars": _stars(pa), "z": float(z)})
    return out


def run_group_comparison(df: pd.DataFrame, group_col: str, value_col: str,
                         order: list[str] | None = None,
                         method: str = "auto") -> StatResult:
    """Pick and run the correct group-comparison test, honestly disclosed.

    method="auto" (default) chooses the family from a normality check:
      - all groups look normal -> PARAMETRIC
          2 groups : Welch's t-test (never assumes equal variance)
          3+ groups: one-way ANOVA + Tukey HSD post-hoc
      - any group non-normal / n<3 -> NON-PARAMETRIC
          2 groups : Mann-Whitney U
          3+ groups: Kruskal-Wallis + Dunn's post-hoc (Holm-adjusted)
    method="parametric" or "nonparametric" forces the family.
    A Levene test also flags unequal variances for the parametric path.
    """
    d = df[[group_col, value_col]].dropna()
    levels = order or list(pd.unique(d[group_col]))
    d = d[d[group_col].isin(levels)]
    samples = [d.loc[d[group_col] == g, value_col].values for g in levels]

    all_normal, warnings = _check_normality(levels, samples)

    if method == "auto":
        use_nonparam = not all_normal
    elif method == "nonparametric":
        use_nonparam = True
    elif method == "parametric":
        use_nonparam = False
    else:
        raise ValueError(f"method must be auto/parametric/nonparametric, got {method!r}")

    if method == "auto" and use_nonparam:
        warnings.append("Data not confirmed normal — used a non-parametric test "
                        "(rank-based). Choose 'parametric' to override.")

    # ---- non-parametric family ------------------------------------------------
    if use_nonparam:
        if len(levels) == 2:
            U, p = stats.mannwhitneyu(samples[0], samples[1], alternative="two-sided")
            return StatResult(
                test_name="Mann-Whitney U",
                summary=f"Mann-Whitney U={U:.0f}, p {fmt_p(p)}",
                omnibus_p=float(p), method="nonparametric",
                posthoc_name="Mann-Whitney U",
                pairwise=[{"g1": levels[0], "g2": levels[1], "p": float(p),
                           "p_text": fmt_p(p), "stars": _stars(p)}],
                warnings=warnings,
            )
        H, p = stats.kruskal(*samples)
        res = StatResult(
            test_name="Kruskal-Wallis",
            summary=f"Kruskal-Wallis H({len(levels)-1})={H:.2f}, p {fmt_p(p)}",
            omnibus_p=float(p), method="nonparametric",
            posthoc_name="Dunn's test (Holm)", warnings=warnings,
        )
        res.pairwise = _dunn_test(d[value_col].to_numpy(dtype=float),
                                  d[group_col].to_numpy(), levels)
        return res

    # ---- parametric family ----------------------------------------------------
    # Homogeneity-of-variance flag (informational; Welch paths tolerate it)
    if all(len(s) >= 2 for s in samples):
        try:
            if stats.levene(*samples, center="median").pvalue < 0.05:
                warnings.append("Group variances differ (Levene p<0.05); for 3+ groups "
                                "a Welch ANOVA / Games-Howell post-hoc would be safer.")
        except Exception:
            pass

    if len(levels) == 2:
        t, p = stats.ttest_ind(samples[0], samples[1], equal_var=False)
        return StatResult(
            test_name="Welch's t-test",
            summary=f"Welch's t-test: t={t:.2f}, p {fmt_p(p)}",
            omnibus_p=float(p), method="parametric", posthoc_name="Welch's t-test",
            pairwise=[{"g1": levels[0], "g2": levels[1], "p": float(p),
                       "p_text": fmt_p(p), "stars": _stars(p)}],
            warnings=warnings,
        )

    F, p = stats.f_oneway(*samples)
    res = StatResult(
        test_name="One-way ANOVA",
        summary=f"one-way ANOVA: F({len(levels)-1},{len(d)-len(levels)})={F:.2f}, p {fmt_p(p)}",
        omnibus_p=float(p), method="parametric",
        posthoc_name="Tukey HSD" if _HAS_SM else "", warnings=warnings,
    )
    if _HAS_SM:
        tuk = pairwise_tukeyhsd(d[value_col], d[group_col])
        tbl = pd.DataFrame(tuk._results_table.data[1:], columns=tuk._results_table.data[0])
        present = set(levels)
        for _, r in tbl.iterrows():
            if str(r["group1"]) not in present or str(r["group2"]) not in present:
                continue
            pv = float(r["p-adj"])
            res.pairwise.append({"g1": str(r["group1"]), "g2": str(r["group2"]),
                                 "p": pv, "p_text": fmt_p(pv), "stars": _stars(pv)})
    return res


# ---------------------------------------------------------------- figure
def make_group_figure(df: pd.DataFrame, group_col: str, value_col: str,
                      order: list[str] | None = None,
                      stat: StatResult | None = None,
                      ylabel: str | None = None, title: str = "",
                      preset: str | None = None, width_mm: float | None = None) -> plt.Figure:
    """Bar (mean±SEM) + individual points + significance brackets. Returns a Figure."""
    apply_preset(preset) if preset else apply_style()
    d = df[[group_col, value_col]].dropna()
    levels = order or list(pd.unique(d[group_col]))
    colors = {g: OKABE_ITO[i % len(OKABE_ITO)] for i, g in enumerate(levels)}

    fig, ax = plt.subplots(figsize=(0.9 * len(levels) + 1.2, 3.0))
    if width_mm:
        set_width_mm(fig, width_mm)
    means, data_max, data_min = [], -np.inf, np.inf
    for i, g in enumerate(levels):
        v = d.loc[d[group_col] == g, value_col].values
        m = float(np.mean(v)); sem = float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else 0.0
        means.append(m)
        data_max = max(data_max, float(np.max(v)), m + sem)
        data_min = min(data_min, float(np.min(v)), m - sem)
        ax.bar(i, m, width=0.6, color=colors[g], alpha=0.55, edgecolor=colors[g], linewidth=1.2, zorder=1)
        ax.errorbar(i, m, yerr=sem, color="black", lw=1.0, capsize=3, zorder=3)
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.28
        ax.scatter(i + jit, v, s=14, facecolor="white", edgecolor=colors[g], linewidth=1.0, zorder=4)

    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, rotation=0)
    ax.set_ylabel(ylabel or value_col.replace("_", " "))
    if title:
        ax.set_title(title, loc="left", fontweight="bold")

    _annotate_brackets(ax, stat, {g: i for i, g in enumerate(levels)}, data_max, data_min)

    if stat:
        fig.text(0.5, -0.04, stat.summary + "   (bars mean±SEM; points = observations)",
                 ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig


def _annotate_brackets(ax, stat: "StatResult | None", idx: dict,
                       data_max: float, data_min: float) -> None:
    """Draw significance brackets above the data, spaced to avoid overlap.

    Placement is driven by the actual data extent (points + error bars), not by
    the bar means — so brackets never collide with tall individual points.
    Baseline is 0 for all-positive data, else the data minimum (handles ΔF/F,
    z-scores, and other signed measures correctly).
    """
    if not (stat and stat.pairwise) or not np.isfinite(data_max):
        # still fix a sensible y-floor for all-positive data
        if np.isfinite(data_min):
            lo = 0.0 if data_min >= 0 else data_min - 0.08 * (data_max - data_min + 1e-9)
            ax.set_ylim(lo, ax.get_ylim()[1])
        return
    span = (data_max - data_min) or (abs(data_max) or 1.0)
    step = 0.09 * span
    top = data_max + 0.10 * span
    drawn = 0
    for pr in stat.pairwise:
        if pr["stars"] == "ns" or pr["g1"] not in idx or pr["g2"] not in idx:
            continue
        x1, x2 = sorted((idx[pr["g1"]], idx[pr["g2"]]))
        y = top + drawn * step
        ax.plot([x1, x1, x2, x2], [y, y + step * 0.3, y + step * 0.3, y], lw=0.9, color="black")
        ax.text((x1 + x2) / 2, y + step * 0.3, pr["stars"], ha="center", va="bottom", fontsize=9)
        drawn += 1
    lo = 0.0 if data_min >= 0 else data_min - 0.08 * span
    ax.set_ylim(lo, top + (drawn + 1) * step)


def make_box_figure(df: pd.DataFrame, group_col: str, value_col: str,
                    order: list[str] | None = None, kind: str = "box",
                    stat: StatResult | None = None,
                    ylabel: str | None = None, title: str = "",
                    preset: str | None = None, width_mm: float | None = None) -> plt.Figure:
    """Box or violin plot with overlaid individual points + significance brackets.

    kind="box"    -> median/IQR box (robust; the honest choice for non-normal data).
    kind="violin" -> kernel-density violin with the same points overlaid.
    Pairs naturally with the non-parametric tests.
    """
    apply_preset(preset) if preset else apply_style()
    d = df[[group_col, value_col]].dropna()
    levels = order or list(pd.unique(d[group_col]))
    d = d[d[group_col].isin(levels)]
    colors = {g: OKABE_ITO[i % len(OKABE_ITO)] for i, g in enumerate(levels)}
    data = [d.loc[d[group_col] == g, value_col].values for g in levels]

    fig, ax = plt.subplots(figsize=(0.9 * len(levels) + 1.2, 3.0))
    if width_mm:
        set_width_mm(fig, width_mm)
    positions = range(len(levels))

    if kind == "violin":
        vp = ax.violinplot(data, positions=list(positions), showextrema=False)
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(colors[levels[i]]); body.set_alpha(0.35)
            body.set_edgecolor(colors[levels[i]]); body.set_linewidth(1.0)
    else:
        bp = ax.boxplot(data, positions=list(positions), widths=0.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", lw=1.2))
        for i, box in enumerate(bp["boxes"]):
            box.set_facecolor(colors[levels[i]]); box.set_alpha(0.35)
            box.set_edgecolor(colors[levels[i]])

    data_max, data_min = -np.inf, np.inf
    for i, g in enumerate(levels):
        v = data[i]
        data_max = max(data_max, float(np.max(v))); data_min = min(data_min, float(np.min(v)))
        jit = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.24
        ax.scatter(i + jit, v, s=14, facecolor="white", edgecolor=colors[g],
                   linewidth=1.0, zorder=4)

    ax.set_xticks(list(positions))
    ax.set_xticklabels(levels, rotation=0)
    ax.set_ylabel(ylabel or value_col.replace("_", " "))
    if title:
        ax.set_title(title, loc="left", fontweight="bold")

    _annotate_brackets(ax, stat, {g: i for i, g in enumerate(levels)}, data_max, data_min)

    if stat:
        cap = "box: median/IQR; points = observations" if kind == "box" \
              else "violin: density; points = observations"
        fig.text(0.5, -0.04, stat.summary + f"   ({cap})",
                 ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- paired data
def run_paired_comparison(before: np.ndarray, after: np.ndarray,
                          method: str = "auto") -> StatResult:
    """Two paired/repeated measures on the same subjects (e.g. pre vs post).

    method="auto": paired t-test if the differences look normal, else Wilcoxon
    signed-rank. Rows with a missing value in either condition are dropped
    pairwise. Using an unpaired test on paired data is a classic error this
    guards against.
    """
    a = np.asarray(before, dtype=float); b = np.asarray(after, dtype=float)
    mask = ~(np.isnan(a) | np.isnan(b))
    a, b = a[mask], b[mask]
    diff = b - a
    warnings: list[str] = []
    if len(diff) < 3:
        warnings.append(f"Only n={len(diff)} complete pairs — results unreliable.")

    normal = True
    if len(diff) >= 3:
        try:
            normal = stats.shapiro(diff).pvalue >= 0.05
        except Exception:
            normal = True
    use_nonparam = (method == "nonparametric") or (method == "auto" and not normal)
    if method not in ("auto", "parametric", "nonparametric"):
        raise ValueError(f"method must be auto/parametric/nonparametric, got {method!r}")

    if use_nonparam:
        try:
            W, p = stats.wilcoxon(a, b)
        except ValueError:  # e.g. all differences zero
            W, p = float("nan"), 1.0
        name = "Wilcoxon signed-rank"
        summary = f"Wilcoxon signed-rank: W={W:.0f}, p {fmt_p(p)}" if W == W \
                  else f"Wilcoxon signed-rank: p {fmt_p(p)}"
        method_used = "nonparametric"
    else:
        t, p = stats.ttest_rel(a, b)
        name = "Paired t-test"
        summary = f"Paired t-test: t={t:.2f}, p {fmt_p(p)}"
        method_used = "parametric"

    return StatResult(
        test_name=name, summary=summary, omnibus_p=float(p),
        method=method_used, posthoc_name=name,
        pairwise=[{"g1": "before", "g2": "after", "p": float(p),
                   "p_text": fmt_p(p), "stars": _stars(p)}],
        warnings=warnings,
    )


def make_paired_figure(before: np.ndarray, after: np.ndarray,
                       labels: tuple[str, str] = ("Before", "After"),
                       stat: StatResult | None = None,
                       ylabel: str = "value", title: str = "",
                       preset: str | None = None, width_mm: float | None = None) -> plt.Figure:
    """Before/after plot: each subject is a faint connecting line; means in bold.

    The connecting lines are the whole point of a paired figure — they show the
    within-subject change a bar chart hides.
    """
    apply_preset(preset) if preset else apply_style()
    a = np.asarray(before, dtype=float); b = np.asarray(after, dtype=float)
    mask = ~(np.isnan(a) | np.isnan(b))
    a, b = a[mask], b[mask]
    c0, c1 = OKABE_ITO[0], OKABE_ITO[1]

    fig, ax = plt.subplots(figsize=(2.6, 3.0))
    if width_mm:
        set_width_mm(fig, width_mm)
    for ai, bi in zip(a, b):
        ax.plot([0, 1], [ai, bi], color="0.6", lw=0.7, alpha=0.7, zorder=1)
    ax.scatter(np.zeros_like(a), a, s=16, facecolor="white", edgecolor=c0, linewidth=1.0, zorder=3)
    ax.scatter(np.ones_like(b), b, s=16, facecolor="white", edgecolor=c1, linewidth=1.0, zorder=3)

    for x, vals, col in ((0, a, c0), (1, b, c1)):
        m = float(np.mean(vals))
        sem = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
        ax.errorbar(x + (0.12 if x else -0.12), m, yerr=sem, fmt="o", color=col,
                    capsize=3, ms=5, lw=1.2, zorder=4)

    ax.set_xticks([0, 1]); ax.set_xticklabels(list(labels))
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")

    if stat and stat.pairwise and stat.pairwise[0]["stars"] != "ns":
        both = np.concatenate([a, b]); hi = float(np.max(both)); lo = float(np.min(both))
        span = (hi - lo) or (abs(hi) or 1.0)
        y = hi + 0.10 * span
        ax.plot([0, 0, 1, 1], [y, y + 0.03 * span, y + 0.03 * span, y], lw=0.9, color="black")
        ax.text(0.5, y + 0.03 * span, stat.pairwise[0]["stars"], ha="center", va="bottom", fontsize=9)
        ax.set_ylim(min(lo, 0) if lo >= 0 else lo - 0.08 * span, y + 0.14 * span)

    if stat:
        fig.text(0.5, -0.04, stat.summary + "   (lines = individual subjects; markers mean±SEM)",
                 ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------- time-course / ΔF/F
def timecourse_from_wide(df: pd.DataFrame, time_col: str,
                         trace_cols: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Wide time-course table -> (time vector, traces matrix of shape n_subjects×n_time).

    The near-universal export shape from calcium/photometry pipelines is one row
    per time point and one column per subject/trial. We transpose to subjects×time
    so every downstream op (mean, SEM, AUC) is 'across subjects at each time'.
    """
    t = np.asarray(df[time_col], dtype=float)
    m = np.asarray(df[trace_cols], dtype=float).T
    return t, m


def timecourse_stats(traces: np.ndarray, error: str = "sem") -> tuple[np.ndarray, np.ndarray]:
    """Across-subject mean and error band at each time point (NaN-aware).

    error='sem' (default) or 'sd'. Uses ddof=1 (sample) and, for SEM, divides by
    sqrt of the *non-NaN* count per column so ragged traces stay correct.
    """
    traces = np.asarray(traces, dtype=float)
    mean = np.nanmean(traces, axis=0)
    sd = np.nanstd(traces, axis=0, ddof=1)
    if error == "sd":
        return mean, sd
    n = np.sum(~np.isnan(traces), axis=0)
    return mean, sd / np.sqrt(np.maximum(n, 1))


def timecourse_auc(time: np.ndarray, traces: np.ndarray,
                   window: tuple[float, float]) -> np.ndarray:
    """Trapezoidal area under each subject's trace within [t0, t1].

    This is the honest way to bring statistics to a time-course: reduce each
    subject to one number in a pre-registered window, then run a standard test
    (run_group_comparison / run_paired_comparison) on those AUCs. Avoids the
    invalid 'test every time-bin' approach.
    """
    time = np.asarray(time, dtype=float)
    traces = np.asarray(traces, dtype=float)
    t0, t1 = window
    mask = (time >= t0) & (time <= t1)
    if mask.sum() < 2:
        raise ValueError("AUC window must cover at least 2 time points.")
    integ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # renamed in NumPy 2.0
    return integ(traces[:, mask], time[mask], axis=1)


def make_timecourse_figure(time: np.ndarray,
                           traces_by_condition: dict,
                           ylabel: str = "ΔF/F (%)", xlabel: str = "Time (s)",
                           error: str = "sem",
                           event_window=None, event_label: str | None = None,
                           title: str = "",
                           preset: str | None = None, width_mm: float | None = None) -> plt.Figure:
    """Mean±SEM (or SD) trace with shaded band; one line per condition.

    traces_by_condition : {label -> array (n_subjects × n_time)}.
    event_window        : a float draws a dashed line at that time (stimulus onset);
                          a (t0, t1) tuple shades the stimulus epoch.
    A baseline at ΔF/F = 0 is drawn because ΔF/F is defined relative to baseline.
    """
    apply_preset(preset) if preset else apply_style()
    time = np.asarray(time, dtype=float)
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    if width_mm:
        set_width_mm(fig, width_mm)

    ax.axhline(0, color="0.7", lw=0.6, zorder=1)
    if event_window is not None:
        if isinstance(event_window, (int, float)):
            ax.axvline(float(event_window), color="0.4", lw=0.8, ls="--", zorder=1)
        else:
            ax.axvspan(float(event_window[0]), float(event_window[1]),
                       color="0.85", zorder=0, label=event_label)

    for i, (label, traces) in enumerate(traces_by_condition.items()):
        mean, err = timecourse_stats(traces, error=error)
        c = OKABE_ITO[i % len(OKABE_ITO)]
        ax.fill_between(time, mean - err, mean + err, color=c, alpha=0.25,
                        linewidth=0, zorder=2)
        ax.plot(time, mean, color=c, lw=1.3, zorder=3, label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(float(time.min()), float(time.max()))
    if title:
        ax.set_title(title, loc="left", fontweight="bold")
    if len(traces_by_condition) > 1 or event_label:
        ax.legend(frameon=False, fontsize=6.5, loc="best")

    fig.text(0.5, -0.04, f"line = mean; shading = ±{error.upper()} across subjects",
             ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig


def figure_to_bytes(fig: plt.Figure, fmt: str = "png", exact: bool = False) -> bytes:
    """Serialize a figure to bytes.

    exact=False (default): trim surrounding whitespace (bbox_inches='tight') —
        best for previews and general use.
    exact=True: preserve the figure's exact canvas size, so a width set via
        width_mm is honoured to the millimetre in the output file. Use this for
        journal submission where the column width must be exact.
    """
    buf = io.BytesIO()
    if exact:
        fig.savefig(buf, format=fmt, facecolor="white")
    else:
        fig.savefig(buf, format=fmt, bbox_inches="tight", facecolor="white")
    buf.seek(0)
    return buf.read()
