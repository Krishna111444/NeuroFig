"""
journal_figures.py — the 20 recurring top-journal figure types, on the shared
NeuroFig engine (Okabe-Ito palette, journal presets, exact-mm export, editing).

Each figure is a thin function that builds a matplotlib Figure from data; a sample
data generator is bundled so every one is demoable and testable with no upload.
The FIGURES registry at the bottom lets the app enumerate them.

Design point: the expensive parts (style, colour, export) live in neurofig_core,
so adding a figure is cheap — the platform is the moat, not any single plot.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle
from matplotlib.path import Path
import matplotlib.patches as mpatches
from scipy import stats

from neurofig_core import OKABE_ITO, apply_preset, apply_style, set_width_mm


def _new(preset, width_mm, w=4.0, h=3.0):
    apply_preset(preset) if preset else apply_style()
    fig, ax = plt.subplots(figsize=(w, h))
    if width_mm:
        set_width_mm(fig, width_mm)
    return fig, ax


def _rng(seed=0):
    return np.random.default_rng(seed)


# ============================================================ 01 individual error dots
def error_dotplot_individual(preset=None, width_mm=None):
    """Per-item estimate ± CI (a forest-style caterpillar of individuals)."""
    fig, ax = _new(preset, width_mm, 3.2, 3.4)
    rng = _rng(1)
    n = 20
    est = np.sort(rng.normal(0, 1, n))
    err = rng.uniform(0.2, 0.5, n)
    y = np.arange(n)
    ax.errorbar(est, y, xerr=err, fmt="o", ms=4, color=OKABE_ITO[0],
                ecolor="0.6", lw=0.8, mfc="white")
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set_yticks([]); ax.set_xlabel("Effect size (± 95% CI)"); ax.set_ylabel("Individuals")
    fig.tight_layout(); return fig


# ============================================================ 02 grouped error dots
def error_dotplot_grouped(preset=None, width_mm=None):
    """Mean ± CI dots for several measures across groups (dodged)."""
    fig, ax = _new(preset, width_mm, 4.2, 3.0)
    rng = _rng(2)
    measures = ["M1", "M2", "M3", "M4"]
    groups = ["Control", "Treated"]
    for gi, g in enumerate(groups):
        means = rng.normal(gi * 0.6, 0.5, len(measures))
        err = rng.uniform(0.15, 0.35, len(measures))
        y = np.arange(len(measures)) + (gi - 0.5) * 0.2
        ax.errorbar(means, y, xerr=err, fmt="o", ms=5, color=OKABE_ITO[gi],
                    ecolor=OKABE_ITO[gi], lw=1.0, mfc="white", label=g)
    ax.axvline(0, color="0.4", ls="--", lw=0.8)
    ax.set_yticks(range(len(measures))); ax.set_yticklabels(measures)
    ax.set_xlabel("Estimate ± CI"); ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); return fig


# ============================================================ 03 volcano
def volcano(preset=None, width_mm=None):
    """Differential-expression volcano: log2 fold-change vs −log10 p."""
    fig, ax = _new(preset, width_mm, 3.6, 3.2)
    rng = _rng(3)
    n = 2000
    lfc = rng.normal(0, 1.2, n)
    p = 10 ** (-np.abs(rng.normal(0, 1.5, n)) - 0.1 * lfc ** 2)
    neglogp = -np.log10(p)
    up = (lfc > 1) & (neglogp > 2)
    dn = (lfc < -1) & (neglogp > 2)
    ns = ~(up | dn)
    ax.scatter(lfc[ns], neglogp[ns], s=4, color="0.7", alpha=0.5, linewidth=0)
    ax.scatter(lfc[up], neglogp[up], s=6, color=OKABE_ITO[1], linewidth=0, label="Up")
    ax.scatter(lfc[dn], neglogp[dn], s=6, color=OKABE_ITO[0], linewidth=0, label="Down")
    ax.axhline(2, color="0.4", ls="--", lw=0.7); ax.axvline(1, color="0.4", ls="--", lw=0.7)
    ax.axvline(-1, color="0.4", ls="--", lw=0.7)
    ax.set_xlabel("log₂ fold change".replace("₂", "2")); ax.set_ylabel("-log10 p")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout(); return fig


# ============================================================ 04 Manhattan
def manhattan(preset=None, width_mm=None):
    """GWAS/TWAS Manhattan: −log10 p along the genome, chromosomes alternating."""
    fig, ax = _new(preset, width_mm, 5.2, 2.6)
    rng = _rng(4)
    chroms = list(range(1, 13))
    sizes = rng.integers(200, 600, len(chroms))
    offset = 0; xt = []
    for ci, (c, s) in enumerate(zip(chroms, sizes)):
        pos = offset + np.arange(s)
        y = -np.log10(rng.uniform(0, 1, s))
        # inject a couple of hits
        hits = rng.choice(s, 3, replace=False)
        y[hits] += rng.uniform(3, 6, 3)
        ax.scatter(pos, y, s=3, color=OKABE_ITO[ci % 2], linewidth=0)
        xt.append(offset + s / 2); offset += s
    ax.axhline(-np.log10(5e-8), color=OKABE_ITO[1], ls="--", lw=0.8)
    ax.set_xticks(xt); ax.set_xticklabels(chroms, fontsize=5)
    ax.set_xlabel("Chromosome"); ax.set_ylabel("-log10 p")
    fig.tight_layout(); return fig


# ============================================================ 05 paired boxplot
def paired_boxplot(preset=None, width_mm=None):
    """Paired boxplots with within-subject connecting lines."""
    fig, ax = _new(preset, width_mm, 2.8, 3.2)
    rng = _rng(5)
    pre = rng.normal(10, 2, 18); post = pre + rng.normal(2.2, 1.0, 18)
    data = [pre, post]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black"))
    for i, box in enumerate(bp["boxes"]):
        box.set_facecolor(OKABE_ITO[i]); box.set_alpha(0.35); box.set_edgecolor(OKABE_ITO[i])
    for a, b in zip(pre, post):
        ax.plot([0, 1], [a, b], color="0.6", lw=0.5, alpha=0.6, zorder=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre", "Post"]); ax.set_ylabel("Value")
    fig.tight_layout(); return fig


# ============================================================ 06 raincloud
def _half_violin(ax, data, pos, color, side=1, width=0.35):
    kde = stats.gaussian_kde(data)
    ys = np.linspace(min(data), max(data), 100)
    dens = kde(ys); dens = dens / dens.max() * width
    ax.fill_betweenx(ys, pos, pos + side * dens, color=color, alpha=0.35, linewidth=0)


def raincloud(preset=None, width_mm=None):
    """Raincloud: half-violin density + boxplot + jittered raw points, per group."""
    fig, ax = _new(preset, width_mm, 4.0, 3.0)
    rng = _rng(6)
    groups = ["A", "B", "C"]
    for i, g in enumerate(groups):
        d = rng.normal(i * 1.2, 1.0, 60)
        c = OKABE_ITO[i]
        _half_violin(ax, d, i + 0.05, c, side=1)
        ax.boxplot(d, positions=[i - 0.12], widths=0.12, patch_artist=True, showfliers=False,
                   medianprops=dict(color="black"),
                   boxprops=dict(facecolor=c, alpha=0.4, edgecolor=c))
        jit = i - 0.28 + (rng.random(len(d)) - 0.5) * 0.12
        ax.scatter(jit, d, s=6, color=c, alpha=0.6, linewidth=0)
    ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups); ax.set_ylabel("Value")
    fig.tight_layout(); return fig


# ============================================================ 09 grouped bars + letters
def _compact_letters(levels, sig_pairs):
    """Very small compact-letter-display: groups sharing a letter are NOT different."""
    letters = {g: set() for g in levels}
    cur = 0
    remaining = set(levels)
    # greedy: build maximal cliques of non-significantly-different groups
    def diff(a, b):
        return (a, b) in sig_pairs or (b, a) in sig_pairs
    assigned = {g: "" for g in levels}
    letter_ord = ord("a")
    groups = list(levels)
    used = []
    for g in groups:
        placed = False
        for grp in used:
            if all(not diff(g, h) for h in grp):
                grp.append(g); placed = True; break
        if not placed:
            used.append([g])
    for i, grp in enumerate(used):
        for g in grp:
            assigned[g] += chr(letter_ord + i)
    return assigned


def grouped_bars_letters(preset=None, width_mm=None):
    """Grouped mean±SEM bars annotated with post-hoc significance letters."""
    fig, ax = _new(preset, width_mm, 3.6, 3.0)
    rng = _rng(9)
    levels = ["Ctrl", "LowD", "HighD", "Combo"]
    data = {g: rng.normal(m, 1.0, 12) for g, m in zip(levels, [5, 6, 9, 9.3])}
    # pairwise Welch tests -> significant pairs
    sig = set()
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            p = stats.ttest_ind(data[levels[i]], data[levels[j]], equal_var=False).pvalue
            if p < 0.05:
                sig.add((levels[i], levels[j]))
    letters = _compact_letters(levels, sig)
    for i, g in enumerate(levels):
        m = data[g].mean(); sem = data[g].std(ddof=1) / np.sqrt(len(data[g]))
        ax.bar(i, m, 0.6, color=OKABE_ITO[i], alpha=0.6, edgecolor=OKABE_ITO[i])
        ax.errorbar(i, m, yerr=sem, color="black", capsize=3, lw=1)
        ax.text(i, m + sem + 0.3, letters[g], ha="center", fontsize=9, fontweight="bold")
    ax.set_xticks(range(len(levels))); ax.set_xticklabels(levels); ax.set_ylabel("Response")
    fig.tight_layout(); return fig


# ============================================================ 13 binned heatmap
def binned_heatmap(preset=None, width_mm=None):
    """Discrete (binned) heatmap — values mapped to a few ordinal colour bins."""
    fig, ax = _new(preset, width_mm, 3.6, 3.0)
    rng = _rng(13)
    M = rng.normal(0, 1, (8, 10))
    from matplotlib.colors import BoundaryNorm, ListedColormap
    bounds = [-3, -1.5, -0.5, 0.5, 1.5, 3]
    cmap = ListedColormap([OKABE_ITO[0], "#9ecae1", "#f0f0f0", "#fdae6b", OKABE_ITO[1]])
    norm = BoundaryNorm(bounds, cmap.N)
    im = ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xlabel("Samples"); ax.set_ylabel("Features")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, ticks=bounds)
    fig.tight_layout(); return fig


# ============================================================ 15 correlation heatmap
def correlation_heatmap(preset=None, width_mm=None):
    """Lower-triangle correlation matrix with a diverging scale."""
    fig, ax = _new(preset, width_mm, 3.4, 3.0)
    rng = _rng(15)
    k = 8
    base = rng.normal(0, 1, (200, k))
    base[:, 1] += base[:, 0]; base[:, 3] -= base[:, 2]  # induce correlations
    C = np.corrcoef(base, rowvar=False)
    mask = np.triu(np.ones_like(C, dtype=bool), k=1)
    Cm = np.ma.array(C, mask=mask)
    im = ax.imshow(Cm, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(k)); ax.set_yticks(range(k))
    ax.set_xticklabels([f"g{i}" for i in range(k)], fontsize=5)
    ax.set_yticklabels([f"g{i}" for i in range(k)], fontsize=5)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label="Pearson r")
    fig.tight_layout(); return fig


# ============================================================ 16 polar heatmap
def polar_heatmap(preset=None, width_mm=None):
    """Circular heatmap: sectors × rings (e.g. samples × features around a circle)."""
    apply_preset(preset) if preset else apply_style()
    fig, ax = plt.subplots(figsize=(3.4, 3.2), subplot_kw=dict(projection="polar"))
    if width_mm:
        set_width_mm(fig, width_mm)
    rng = _rng(16)
    n_sec, n_ring = 24, 5
    M = rng.normal(0, 1, (n_ring, n_sec))
    theta = np.linspace(0, 2 * np.pi, n_sec + 1)
    r = np.linspace(1, 2, n_ring + 1)
    T, R = np.meshgrid(theta, r)
    ax.pcolormesh(T, R, M, cmap="viridis", shading="flat")
    ax.set_yticklabels([]); ax.set_xticklabels([]); ax.spines["polar"].set_visible(False)
    fig.tight_layout(); return fig


# ============================================================ 20 split violin
def split_violin(preset=None, width_mm=None):
    """Split violin: two subgroups compared as mirrored half-violins per category."""
    fig, ax = _new(preset, width_mm, 4.0, 3.0)
    rng = _rng(20)
    cats = ["Region 1", "Region 2", "Region 3"]
    for i, cat in enumerate(cats):
        a = rng.normal(i, 1.0, 80); b = rng.normal(i + 0.6, 1.1, 80)
        _half_violin(ax, a, i, OKABE_ITO[0], side=-1, width=0.4)
        _half_violin(ax, b, i, OKABE_ITO[1], side=1, width=0.4)
    ax.set_xticks(range(len(cats))); ax.set_xticklabels(cats); ax.set_ylabel("Value")
    ax.legend(handles=[mpatches.Patch(color=OKABE_ITO[0], alpha=0.4, label="Sham"),
                       mpatches.Patch(color=OKABE_ITO[1], alpha=0.4, label="Lesion")],
              frameon=False, fontsize=6.5)
    fig.tight_layout(); return fig


# ============================================================ 07 swimmer plot
def swimmer_plot(preset=None, width_mm=None):
    """Swimmer plot: one horizontal bar per subject (time on study), with a
    response marker and an arrow for ongoing subjects (oncology staple)."""
    fig, ax = _new(preset, width_mm, 3.6, 3.4)
    rng = _rng(7)
    n = 16
    dur = np.sort(rng.uniform(2, 24, n))
    resp = rng.integers(0, 3, n)
    ongoing = rng.random(n) > 0.5
    cols = [OKABE_ITO[0], OKABE_ITO[4], OKABE_ITO[2]]
    for i in range(n):
        ax.barh(i, dur[i], height=0.6, color=cols[resp[i]], alpha=0.8)
        ax.scatter(dur[i] * rng.uniform(0.2, 0.5), i, marker="o", s=14, color="black", zorder=3)
        if ongoing[i]:
            ax.annotate("", xy=(dur[i] + 1.2, i), xytext=(dur[i], i),
                        arrowprops=dict(arrowstyle="->", lw=1.0))
    ax.set_yticks([]); ax.set_xlabel("Months on study"); ax.set_ylabel("Patients")
    ax.legend(handles=[mpatches.Patch(color=cols[k], label=l)
                       for k, l in enumerate(["PD", "SD", "PR/CR"])],
              frameon=False, fontsize=6, loc="lower right")
    fig.tight_layout(); return fig


# ============================================================ 08 importance + streams
def importance_streams(preset=None, width_mm=None):
    """Feature-importance bars beside a streamgraph of those features over time."""
    apply_preset(preset) if preset else apply_style()
    fig, (axl, axr) = plt.subplots(1, 2, figsize=(5.4, 3.0),
                                   gridspec_kw=dict(width_ratios=[1, 1.4]))
    if width_mm:
        set_width_mm(fig, width_mm)
    rng = _rng(8)
    feats = [f"f{i}" for i in range(7)]
    imp = np.sort(rng.uniform(0.05, 1, len(feats)))
    axl.barh(range(len(feats)), imp, color=OKABE_ITO[0], alpha=0.8)
    axl.set_yticks(range(len(feats))); axl.set_yticklabels(feats, fontsize=6)
    axl.set_xlabel("Importance"); axl.invert_yaxis()
    t = np.linspace(0, 10, 60)
    streams = np.abs(rng.normal(1, 0.3, (len(feats), t.size)).cumsum(axis=1))
    streams = streams / streams.sum(0)
    base = -streams.sum(0) / 2
    for i in range(len(feats)):
        axr.fill_between(t, base, base + streams[i], color=OKABE_ITO[i % len(OKABE_ITO)],
                         alpha=0.8, linewidth=0)
        base = base + streams[i]
    axr.set_xlabel("Time"); axr.set_yticks([]); axr.set_ylabel("Relative abundance")
    fig.tight_layout(); return fig


# ---------- geometry helpers for chords / ribbons ----------
def _chord(ax, p0, p1, color, lw=0.8, alpha=0.6):
    verts = [p0, (0, 0), p1]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc="none", ec=color, lw=lw, alpha=alpha))


def _ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, color, alpha=0.55):
    cx = (x0 + x1) / 2
    verts = [(x0, y0a), (cx, y0a), (cx, y1a), (x1, y1a),
             (x1, y1b), (cx, y1b), (cx, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), fc=color, ec="none", alpha=alpha))


# ============================================================ 10 circos fusion-gene
def circos(preset=None, width_mm=None):
    """Circular ideogram of chromosomes with fusion-gene chords across the circle."""
    fig, ax = _new(preset, width_mm, 3.6, 3.6)
    rng = _rng(10)
    n = 10
    sizes = rng.uniform(0.6, 1.4, n); sizes = sizes / sizes.sum() * 360
    start = 90.0; segs = []
    for i in range(n):
        a0, a1 = start, start - sizes[i]
        ax.add_patch(Wedge((0, 0), 1.0, a1, a0, width=0.12,
                     facecolor=OKABE_ITO[i % len(OKABE_ITO)], edgecolor="white", lw=0.5))
        segs.append((a0, a1)); start = a1

    def pt(seg):
        a = np.deg2rad(rng.uniform(min(seg), max(seg)))
        return (0.86 * np.cos(a), 0.86 * np.sin(a))
    for _ in range(8):
        i, j = int(rng.integers(0, n)), int(rng.integers(0, n))
        if i != j:
            _chord(ax, pt(segs[i]), pt(segs[j]), OKABE_ITO[3], lw=0.8, alpha=0.5)
    ax.set_xlim(-1.15, 1.15); ax.set_ylim(-1.15, 1.15); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(); return fig


# ============================================================ 11 module network
def network(preset=None, width_mm=None):
    """Module-interaction network: nodes coloured by module, spring layout."""
    fig, ax = _new(preset, width_mm, 3.6, 3.4)
    rng = _rng(11)
    n = 28; modules = rng.integers(0, 4, n)
    edges = set()
    while len(edges) < 46:
        a, b = int(rng.integers(0, n)), int(rng.integers(0, n))
        if a != b and (modules[a] == modules[b] or rng.random() < 0.3):
            edges.add((min(a, b), max(a, b)))
    edges = list(edges)
    pos = rng.normal(0, 1, (n, 2))
    for _ in range(60):
        disp = np.zeros((n, 2))
        for i in range(n):
            d = pos[i] - pos; dist = np.linalg.norm(d, axis=1); dist[i] = 1
            disp[i] += (d / dist[:, None] ** 2).sum(0) * 0.05
        for a, b in edges:
            d = pos[a] - pos[b]; disp[a] -= d * 0.05; disp[b] += d * 0.05
        pos += np.clip(disp, -0.1, 0.1)
    for a, b in edges:
        ax.plot([pos[a, 0], pos[b, 0]], [pos[a, 1], pos[b, 1]], color="0.75", lw=0.5, zorder=1)
    ax.scatter(pos[:, 0], pos[:, 1], c=[OKABE_ITO[m] for m in modules], s=40,
               edgecolor="white", lw=0.6, zorder=2)
    ax.set_aspect("equal"); ax.axis("off"); fig.tight_layout(); return fig


# ============================================================ 12 sankey + bubble
def sankey_bubble(preset=None, width_mm=None):
    """Two-layer Sankey beside an enrichment bubble panel."""
    apply_preset(preset) if preset else apply_style()
    fig, (axs, axb) = plt.subplots(1, 2, figsize=(5.6, 3.0),
                                   gridspec_kw=dict(width_ratios=[1.3, 1]))
    if width_mm:
        set_width_mm(fig, width_mm)
    rng = _rng(12)
    src = [0.5, 0.3, 0.2]; tgt = [0.4, 0.35, 0.25]
    sy0 = np.concatenate([[0], np.cumsum(src)]); ty0 = np.concatenate([[0], np.cumsum(tgt)])
    flow = rng.uniform(0.02, 0.1, (3, 3))
    so = sy0[:3].copy(); to = ty0[:3].copy()
    for i in range(3):
        for j in range(3):
            f = flow[i, j]
            _ribbon(axs, 0, so[i], so[i] + f, 1, to[j], to[j] + f, OKABE_ITO[i])
            so[i] += f; to[j] += f
    for i in range(3):
        axs.add_patch(Rectangle((-0.03, sy0[i]), 0.03, src[i], color=OKABE_ITO[i]))
        axs.add_patch(Rectangle((1.0, ty0[i]), 0.03, tgt[i], color="0.5"))
    axs.set_xlim(-0.1, 1.1); axs.set_ylim(0, 1.05); axs.axis("off")
    axb.scatter(-np.log10(rng.uniform(1e-6, 0.05, 6)), range(6),
                s=rng.uniform(20, 120, 6), color=OKABE_ITO[1], alpha=0.7)
    axb.set_yticks(range(6)); axb.set_yticklabels([f"pathway {i}" for i in range(6)], fontsize=5)
    axb.set_xlabel("-log10 p"); fig.tight_layout(); return fig


# ============================================================ 14 Mantel composite
def mantel_heatmap(preset=None, width_mm=None):
    """Correlation heatmap with Mantel-style curved links from drivers to the matrix."""
    fig, ax = _new(preset, width_mm, 4.4, 3.2)
    rng = _rng(14)
    k = 6; base = rng.normal(0, 1, (150, k)); base[:, 1] += base[:, 0]
    C = np.corrcoef(base, rowvar=False)
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1, extent=[0, k, 0, k])
    ax.set_xticks(np.arange(k) + 0.5); ax.set_yticks(np.arange(k) + 0.5)
    ax.set_xticklabels([f"t{i}" for i in range(k)], fontsize=5)
    ax.set_yticklabels([f"t{i}" for i in range(k)], fontsize=5)
    for di, d in enumerate(["pH", "Temp", "Depth"]):
        x0, y0 = -1.5, k - 1 - di * 2
        ax.text(x0 - 0.2, y0, d, ha="right", va="center", fontsize=6)
        tgt = int(rng.integers(0, k)); r = rng.uniform(-0.7, 0.7)
        col = OKABE_ITO[1] if r > 0 else OKABE_ITO[0]
        verts = [(x0, y0), (x0 / 2, k / 2), (tgt + 0.5, k)]
        ax.add_patch(mpatches.PathPatch(Path(verts, [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
                     fc="none", ec=col, lw=0.6 + 2 * abs(r), alpha=0.6))
    ax.set_xlim(-3, k); ax.set_ylim(0, k)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="r")
    fig.tight_layout(); return fig


# ============================================================ 17 multi-level sankey
def multilevel_sankey(preset=None, width_mm=None):
    """Three-layer Sankey (hierarchical flow)."""
    fig, ax = _new(preset, width_mm, 5.0, 3.0)
    layers = [[0.5, 0.5], [0.34, 0.33, 0.33], [0.4, 0.35, 0.25]]
    xs = [0, 1, 2]
    ypos = [np.concatenate([[0], np.cumsum(l)]) for l in layers]
    for li in range(len(layers) - 1):
        a = ypos[li][:len(layers[li])].copy(); b = ypos[li + 1][:len(layers[li + 1])].copy()
        for i in range(len(layers[li])):
            for j in range(len(layers[li + 1])):
                f = layers[li][i] * layers[li + 1][j]
                _ribbon(ax, xs[li], a[i], a[i] + f, xs[li + 1], b[j], b[j] + f,
                        OKABE_ITO[i % len(OKABE_ITO)]); a[i] += f; b[j] += f
    for li, l in enumerate(layers):
        for i in range(len(l)):
            ax.add_patch(Rectangle((xs[li] - 0.02, ypos[li][i]), 0.04, l[i], color="0.4"))
    ax.set_xlim(-0.2, 2.2); ax.set_ylim(0, 1.05); ax.axis("off")
    fig.tight_layout(); return fig


# ============================================================ 18 treemap
def _squarify(sizes, x, y, dx, dy):
    sizes = [s * dx * dy / sum(sizes) for s in sizes]
    sizes = sorted(sizes, reverse=True)
    rects = []

    def worst(row, side):
        s = sum(row)
        return max(max(side * side * r / (s * s) for r in row), s * s / (side * side * min(row)))

    def layout(row, x, y, dx, dy):
        cov = sum(row); rs = []
        if dx >= dy:
            w = cov / dy; yy = y
            for s in row:
                rs.append((x, yy, w, s / w)); yy += s / w
            return rs, x + w, y, dx - w, dy
        h = cov / dx; xx = x
        for s in row:
            rs.append((xx, y, s / h, h)); xx += s / h
        return rs, x, y + h, dx, dy - h

    row = []; i = 0
    while i < len(sizes):
        side = min(dx, dy)
        if not row or worst(row + [sizes[i]], side) <= worst(row, side):
            row.append(sizes[i]); i += 1
        else:
            placed, x, y, dx, dy = layout(row, x, y, dx, dy); rects += placed; row = []
    if row:
        placed, x, y, dx, dy = layout(row, x, y, dx, dy); rects += placed
    return rects


def treemap(preset=None, width_mm=None):
    """Squarified treemap of hierarchical proportions."""
    fig, ax = _new(preset, width_mm, 3.8, 3.0)
    rng = _rng(18)
    vals = sorted(rng.uniform(1, 10, 10), reverse=True)
    for i, (x, y, w, h) in enumerate(_squarify(vals, 0, 0, 1, 1)):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=OKABE_ITO[i % len(OKABE_ITO)],
                     edgecolor="white", lw=1, alpha=0.85))
        if w > 0.08 and h > 0.08:
            ax.text(x + w / 2, y + h / 2, f"c{i}", ha="center", va="center", fontsize=6)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    fig.tight_layout(); return fig


# ============================================================ 19 mosaic + sunburst
def mosaic_sunburst(preset=None, width_mm=None):
    """Mosaic (marimekko) plus a sunburst of the same 2-level hierarchy."""
    apply_preset(preset) if preset else apply_style()
    fig = plt.figure(figsize=(5.6, 3.0))
    if width_mm:
        set_width_mm(fig, width_mm)
    axm = fig.add_subplot(1, 2, 1)
    axs = fig.add_subplot(1, 2, 2, projection="polar")
    rng = _rng(19)
    x0 = 0
    for ci, w in enumerate([0.4, 0.35, 0.25]):
        parts = rng.uniform(0.2, 1, 3); parts /= parts.sum(); y0 = 0
        for pi, p in enumerate(parts):
            axm.add_patch(Rectangle((x0, y0), w, p, facecolor=OKABE_ITO[pi],
                          edgecolor="white", lw=1, alpha=0.85)); y0 += p
        x0 += w
    axm.set_xlim(0, 1); axm.set_ylim(0, 1); axm.set_title("Mosaic", fontsize=7, loc="left")
    axm.set_xticks([]); axm.set_yticks([])
    a0 = 0
    for ci, w in enumerate(np.array([0.5, 0.3, 0.2]) * 2 * np.pi):
        axs.bar(a0 + w / 2, 1.0, width=w, bottom=1, color=OKABE_ITO[ci], alpha=0.85,
                edgecolor="white")
        sub = rng.uniform(0.2, 1, 3); sub = sub / sub.sum() * w; s0 = a0
        for si, sw in enumerate(sub):
            axs.bar(s0 + sw / 2, 1.0, width=sw, bottom=2, color=OKABE_ITO[si], alpha=0.6,
                    edgecolor="white"); s0 += sw
        a0 += w
    axs.set_xticks([]); axs.set_yticks([]); axs.spines["polar"].set_visible(False)
    axs.set_title("Sunburst", fontsize=7, loc="left")
    fig.tight_layout(); return fig


# ============================================================ registry
FIGURES = {
    "01 Individual error dot-plot": error_dotplot_individual,
    "02 Grouped error dot-plot": error_dotplot_grouped,
    "03 Volcano": volcano,
    "04 Manhattan (GWAS/TWAS)": manhattan,
    "05 Paired boxplot": paired_boxplot,
    "06 Raincloud": raincloud,
    "07 Swimmer plot": swimmer_plot,
    "08 Importance + streams": importance_streams,
    "09 Grouped bars + letters": grouped_bars_letters,
    "10 Circos (fusion genes)": circos,
    "11 Module network": network,
    "12 Sankey + enrichment bubbles": sankey_bubble,
    "13 Binned heatmap": binned_heatmap,
    "14 Mantel composite heatmap": mantel_heatmap,
    "15 Correlation heatmap": correlation_heatmap,
    "16 Polar heatmap": polar_heatmap,
    "17 Multi-level Sankey": multilevel_sankey,
    "18 Treemap": treemap,
    "19 Mosaic + sunburst": mosaic_sunburst,
    "20 Split violin": split_violin,
}


# ============================================================================
# REAL-DATA VERSIONS (upload-driven), hardened by inversion: for each figure we
# enumerated how it could break or mislead, and guard every case. Each returns
# (figure, warnings) so the UI can surface what was cleaned or dropped.
# ============================================================================

def _coerce(df, cols):
    out = df[cols].copy()
    for c in cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


# ---------------------------------------------------------------- volcano (real)
def volcano_from_data(df, fc_col, p_col, *, fc_is_log2=True, p_is_neglog10=False,
                      fc_thresh=1.0, p_thresh=0.05, label_col=None, top_n=8,
                      preset=None, width_mm=None):
    """Volcano from a differential-expression table.

    Inversion guards:
      • p = 0            -> floored to the smallest non-zero p (else -log10 = ∞)
      • p outside [0,1]  -> dropped (invalid probabilities)
      • linear FC ≤ 0    -> dropped (log2 undefined)
      • NaN/non-numeric  -> dropped
      • already -log10 p -> p_is_neglog10 avoids a double transform
    """
    warnings = []
    keep_cols = [fc_col, p_col] + ([label_col] if label_col else [])
    d = _coerce(df, [fc_col, p_col]).assign(
        **({label_col: df[label_col]} if label_col else {}))
    n0 = len(d)
    d = d.dropna(subset=[fc_col, p_col])
    if len(d) < n0:
        warnings.append(f"Dropped {n0 - len(d)} rows with missing/non-numeric values.")

    fc = d[fc_col].to_numpy(float)
    if not fc_is_log2:
        bad = fc <= 0
        if bad.any():
            warnings.append(f"Dropped {int(bad.sum())} rows with non-positive linear "
                            f"fold-change (log2 undefined).")
            d = d[~bad]; fc = d[fc_col].to_numpy(float)
        fc = np.log2(fc)

    praw = d[p_col].to_numpy(float)
    if p_is_neglog10:
        neglogp = praw
    else:
        bad = (praw < 0) | (praw > 1)
        if bad.any():
            warnings.append(f"Dropped {int(bad.sum())} p-values outside [0, 1].")
            d = d[~bad]; fc = fc[~bad]; praw = praw[~bad]
        nz = praw[praw > 0]
        floor = nz.min() if nz.size else 1e-300
        nzero = int((praw == 0).sum())
        if nzero:
            warnings.append(f"{nzero} p-values equal 0 were floored to the smallest "
                            f"non-zero p ({floor:.1e}) to avoid infinite -log10.")
        praw = np.where(praw == 0, floor, praw)
        neglogp = -np.log10(praw)

    if len(fc) == 0:
        raise ValueError("No valid rows remain after cleaning — check the columns chosen.")

    p_line = p_thresh if p_is_neglog10 else -np.log10(p_thresh)
    up = (fc > fc_thresh) & (neglogp > p_line)
    dn = (fc < -fc_thresh) & (neglogp > p_line)
    ns = ~(up | dn)

    fig, ax = _new(preset, width_mm, 3.6, 3.4)
    ax.scatter(fc[ns], neglogp[ns], s=5, color="0.7", alpha=0.5, linewidth=0)
    ax.scatter(fc[up], neglogp[up], s=8, color=OKABE_ITO[1], linewidth=0, label="Up")
    ax.scatter(fc[dn], neglogp[dn], s=8, color=OKABE_ITO[0], linewidth=0, label="Down")
    ax.axhline(p_line, color="0.4", ls="--", lw=0.7)
    ax.axvline(fc_thresh, color="0.4", ls="--", lw=0.7)
    ax.axvline(-fc_thresh, color="0.4", ls="--", lw=0.7)
    if label_col and top_n:
        score = neglogp * (up | dn)
        idx = np.argsort(score)[-int(top_n):]
        labs = d[label_col].to_numpy()
        for k in idx:
            if up[k] or dn[k]:
                ax.annotate(str(labs[k]), (fc[k], neglogp[k]), fontsize=5,
                            xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 p")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    fig.tight_layout()
    return fig, warnings


# ---------------------------------------------------------------- raincloud (real)
def raincloud_from_data(df, group_col, value_col, order=None,
                        preset=None, width_mm=None):
    """Raincloud from long data (group column + numeric value column).

    Inversion guards:
      • group with n < 2 or zero spread -> KDE is singular; fall back to points
        (and box when n ≥ 2) instead of crashing
      • non-numeric values -> coerced then dropped
      • empty / single group still renders what it can
    """
    warnings = []
    d = df[[group_col, value_col]].copy()
    d[value_col] = pd.to_numeric(d[value_col], errors="coerce")
    n0 = len(d); d = d.dropna()
    if len(d) < n0:
        warnings.append(f"Dropped {n0 - len(d)} rows with missing/non-numeric values.")
    levels = order or list(dict.fromkeys(d[group_col].tolist()))
    levels = [g for g in levels if (d[group_col] == g).sum() > 0]
    if not levels:
        raise ValueError("No groups with data to plot.")

    fig, ax = _new(preset, width_mm, max(3.0, 0.9 * len(levels) + 1.5), 3.0)
    for i, g in enumerate(levels):
        vals = d.loc[d[group_col] == g, value_col].to_numpy(float)
        c = OKABE_ITO[i % len(OKABE_ITO)]
        can_kde = len(vals) >= 2 and np.std(vals) > 1e-9
        if can_kde:
            try:
                _half_violin(ax, vals, i + 0.05, c, side=1)
                ax.boxplot(vals, positions=[i - 0.12], widths=0.12, patch_artist=True,
                           showfliers=False, medianprops=dict(color="black"),
                           boxprops=dict(facecolor=c, alpha=0.4, edgecolor=c))
            except Exception:
                can_kde = False
        if not can_kde:
            warnings.append(f"Group '{g}' (n={len(vals)}) has too little spread for a "
                            f"density/box — showing points only.")
        rng = np.random.default_rng(i)
        jit = i - 0.28 + (rng.random(len(vals)) - 0.5) * 0.12
        ax.scatter(jit, vals, s=8, color=c, alpha=0.6, linewidth=0)
    ax.set_xticks(range(len(levels))); ax.set_xticklabels(levels)
    ax.set_ylabel(value_col.replace("_", " "))
    fig.tight_layout()
    return fig, warnings


# ---------------------------------------------------------------- correlation (real)
def correlation_heatmap_from_data(df, value_cols, method="pearson", annotate=True,
                                  preset=None, width_mm=None):
    """Correlation heatmap from selected numeric columns.

    Inversion guards:
      • constant column -> correlation is undefined (NaN); dropped with a warning
      • < 2 usable columns -> clear error
      • missing values -> pandas pairwise-complete correlation (no silent NaN wipe)
      • too few complete rows -> warned (estimates unstable)
      • method: pearson (linear) / spearman (monotonic) / kendall
    """
    warnings = []
    if method not in ("pearson", "spearman", "kendall"):
        raise ValueError("method must be pearson, spearman, or kendall.")
    d = _coerce(df, value_cols)
    const = [c for c in value_cols if d[c].nunique(dropna=True) <= 1]
    if const:
        warnings.append(f"Dropped constant column(s) with no variance: {', '.join(map(str, const))}.")
        d = d.drop(columns=const)
    if d.shape[1] < 2:
        raise ValueError("Need at least 2 non-constant numeric columns to correlate.")
    if len(d.dropna()) < 3:
        warnings.append("Fewer than 3 complete rows — correlations may be unstable.")

    C = d.corr(method=method)
    k = C.shape[0]
    fig, ax = _new(preset, width_mm, max(3.0, 0.4 * k + 1.5), max(2.8, 0.4 * k + 1.3))
    im = ax.imshow(C.to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(k)); ax.set_yticks(range(k))
    ax.set_xticklabels(C.columns, rotation=90, fontsize=5)
    ax.set_yticklabels(C.columns, fontsize=5)
    if annotate and k <= 12:
        for i in range(k):
            for j in range(k):
                v = C.iat[i, j]
                if v == v:
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.5,
                            color="white" if abs(v) > 0.6 else "black")
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label=f"{method.title()} r")
    fig.tight_layout()
    return fig, warnings


# ---------------------------------------------------------------- sample templates
def sample_volcano_df():
    rng = np.random.default_rng(3)
    n = 800
    lfc = rng.normal(0, 1.2, n)
    p = 10 ** (-np.abs(rng.normal(0, 1.5, n)) - 0.1 * lfc ** 2)
    p[rng.integers(0, n, 5)] = 0.0  # realistic: some tools report p = 0
    return pd.DataFrame({"gene": [f"G{i}" for i in range(n)],
                         "log2FoldChange": np.round(lfc, 3),
                         "pvalue": p})


def sample_raincloud_df():
    rng = np.random.default_rng(6)
    rows = []
    for g, m in [("A", 0), ("B", 1.2), ("C", 2.4)]:
        for _ in range(60):
            rows.append({"group": g, "value": round(float(rng.normal(m, 1.0)), 3)})
    return pd.DataFrame(rows)


def sample_correlation_df():
    rng = np.random.default_rng(15)
    x = rng.normal(0, 1, (200, 6))
    x[:, 1] += x[:, 0]; x[:, 3] -= x[:, 2]
    return pd.DataFrame(np.round(x, 3), columns=[f"var{i}" for i in range(6)])


REAL_DATA_FIGURES = {
    "03 Volcano": "volcano",
    "06 Raincloud": "raincloud",
    "15 Correlation heatmap": "correlation",
}


# ============================================================================
# MULTIOMICS / SET-COMPARISON FIGURES (Venn, UpSet, enrichment, PCA, composition)
# Each has a bundled sample generator (for the gallery) and, where useful, a
# from_data() for real uploads, plus a sample-CSV maker.
# ============================================================================
from matplotlib.patches import Circle  # noqa: E402


def _set_regions(sets):
    """Region counts for 2 or 3 sets (labelled by membership tuple)."""
    if len(sets) == 2:
        a, b = sets
        return {"10": len(a - b), "01": len(b - a), "11": len(a & b)}
    a, b, c = sets
    return {"100": len(a - b - c), "010": len(b - a - c), "001": len(c - a - b),
            "110": len((a & b) - c), "101": len((a & c) - b), "011": len((b & c) - a),
            "111": len(a & b & c)}


def _draw_venn(ax, sets, labels, colors):
    if len(sets) == 2:
        centres = [(-0.35, 0), (0.35, 0)]; r = 0.6
        pos = {"10": (-0.55, 0), "01": (0.55, 0), "11": (0, 0)}
        lab = {0: (-0.5, 0.65), 1: (0.5, 0.65)}
    else:
        centres = [(-0.35, 0.22), (0.35, 0.22), (0.0, -0.38)]; r = 0.62
        pos = {"100": (-0.55, 0.55), "010": (0.55, 0.55), "001": (0.0, -0.72),
               "110": (0.0, 0.55), "101": (-0.38, -0.18), "011": (0.38, -0.18),
               "111": (0.0, 0.08)}
        lab = {0: (-0.7, 0.7), 1: (0.7, 0.7), 2: (0.0, -0.95)}
    for i, (cx, cy) in enumerate(centres):
        ax.add_patch(Circle((cx, cy), r, facecolor=colors[i], alpha=0.35,
                            edgecolor=colors[i], lw=1.4))
        ax.text(*lab[i], labels[i], color=colors[i], fontweight="bold", ha="center", fontsize=8)
    counts = _set_regions(sets)
    for key, (x, y) in pos.items():
        ax.text(x, y, str(counts.get(key, 0)), ha="center", va="center", fontsize=8)
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.25, 1.15); ax.set_aspect("equal"); ax.axis("off")


def venn(preset=None, width_mm=None):
    """Venn diagram of overlapping feature sets (e.g. genes across omics layers)."""
    fig, ax = _new(preset, width_mm, 3.4, 3.2)
    rng = _rng(30)
    allg = np.array([f"G{i}" for i in range(220)])
    sets = [set(rng.choice(allg, n, replace=False)) for n in (85, 95, 75)]
    _draw_venn(ax, sets, ["RNA-seq", "Proteomics", "ATAC-seq"], OKABE_ITO[:3])
    fig.tight_layout(); return fig


def venn_from_data(df, set_cols, preset=None, width_mm=None):
    """Venn from 2–3 columns, each listing the members of one set."""
    if not (2 <= len(set_cols) <= 3):
        raise ValueError("Venn needs 2 or 3 set columns (use UpSet for more).")
    sets = [set(df[c].dropna().astype(str)) for c in set_cols]
    fig, ax = _new(preset, width_mm, 3.4, 3.2)
    _draw_venn(ax, sets, [str(c) for c in set_cols], OKABE_ITO[:len(sets)])
    fig.tight_layout(); return fig, []


def upset(preset=None, width_mm=None):
    """UpSet plot: intersection sizes across many sets (Venn beyond 3 sets)."""
    apply_preset(preset) if preset else apply_style()
    rng = _rng(31)
    names = ["RNA", "Protein", "ATAC", "Methyl", "Metab"]
    allg = np.arange(400)
    sets = {n: set(rng.choice(allg, s, replace=False)) for n, s in zip(names, [160, 150, 140, 120, 110])}
    import itertools
    combos = []
    for r in range(1, len(names) + 1):
        for combo in itertools.combinations(names, r):
            inter = set(allg)
            for n in combo:
                inter &= sets[n]
            for n in names:
                if n not in combo:
                    inter -= sets[n]
            if len(inter):
                combos.append((combo, len(inter)))
    combos.sort(key=lambda x: -x[1]); combos = combos[:12]

    fig = plt.figure(figsize=(5.2, 3.4))
    if width_mm:
        set_width_mm(fig, width_mm)
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1.2], hspace=0.05)
    axb = fig.add_subplot(gs[0]); axm = fig.add_subplot(gs[1], sharex=axb)
    xs = np.arange(len(combos))
    axb.bar(xs, [c[1] for c in combos], color=OKABE_ITO[0], width=0.6)
    axb.set_ylabel("Intersection"); axb.spines["bottom"].set_visible(False)
    axb.tick_params(labelbottom=False)
    for j, n in enumerate(names):
        for i, (combo, _) in enumerate(combos):
            on = n in combo
            axm.scatter(i, j, s=26, color=OKABE_ITO[0] if on else "0.85", zorder=3)
        # connect the dots of each intersection
    for i, (combo, _) in enumerate(combos):
        ys = [names.index(n) for n in combo]
        if len(ys) > 1:
            axm.plot([i, i], [min(ys), max(ys)], color=OKABE_ITO[0], lw=1.2, zorder=2)
    axm.set_yticks(range(len(names))); axm.set_yticklabels(names, fontsize=6)
    axm.set_xticks([]); axm.set_ylim(-0.6, len(names) - 0.4)
    for sp in ("top", "right", "bottom"):
        axm.spines[sp].set_visible(False)
    fig.tight_layout(); return fig


def enrichment_dotplot(preset=None, width_mm=None):
    """GO/KEGG enrichment dot-plot: gene ratio vs term, size=count, colour=p.adjust."""
    fig, ax = _new(preset, width_mm, 4.2, 3.4)
    rng = _rng(32)
    terms = ["Synaptic signaling", "Ion transport", "Immune response",
             "Cell cycle", "Oxidative phosphorylation", "Apoptosis",
             "MAPK cascade", "Lipid metabolism", "DNA repair", "Chemotaxis"]
    ratio = np.sort(rng.uniform(0.05, 0.45, len(terms)))
    count = rng.integers(8, 80, len(terms))
    padj = 10 ** (-rng.uniform(1.5, 6, len(terms)))
    sc = ax.scatter(ratio, range(len(terms)), s=count * 2.5,
                    c=-np.log10(padj), cmap="viridis", edgecolor="0.3", linewidth=0.4)
    ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=6)
    ax.set_xlabel("Gene ratio")
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02); cb.set_label("-log10 p.adjust", fontsize=6)
    fig.tight_layout(); return fig


def pca_ordination(preset=None, width_mm=None):
    """PCA ordination: samples projected on PC1/PC2, coloured by group."""
    fig, ax = _new(preset, width_mm, 3.8, 3.2)
    rng = _rng(33)
    groups = np.repeat(["Control", "Disease", "Treated"], 12)
    centres = {"Control": 0, "Disease": 4, "Treated": 2}
    X = np.array([rng.normal(centres[g], 1, 40) + rng.normal(0, 1, 40) for g in groups])
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    scores = U * S
    var = S ** 2 / (S ** 2).sum()
    for i, g in enumerate(["Control", "Disease", "Treated"]):
        m = groups == g
        ax.scatter(scores[m, 0], scores[m, 1], s=26, color=OKABE_ITO[i], label=g,
                   edgecolor="white", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({var[0]*100:.0f}%)"); ax.set_ylabel(f"PC2 ({var[1]*100:.0f}%)")
    ax.axhline(0, color="0.85", lw=0.6); ax.axvline(0, color="0.85", lw=0.6)
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout(); return fig


def stacked_composition(preset=None, width_mm=None):
    """Stacked composition bars (relative abundance) — microbiome / cell types."""
    fig, ax = _new(preset, width_mm, 4.2, 3.0)
    rng = _rng(34)
    samples = [f"S{i+1}" for i in range(10)]
    taxa = ["Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria", "Other"]
    comp = rng.uniform(0.2, 1, (len(samples), len(taxa)))
    comp = comp / comp.sum(1, keepdims=True)
    bottom = np.zeros(len(samples))
    for j, tx in enumerate(taxa):
        ax.bar(range(len(samples)), comp[:, j], bottom=bottom,
               color=OKABE_ITO[j % len(OKABE_ITO)], label=tx, width=0.8)
        bottom += comp[:, j]
    ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples, fontsize=6)
    ax.set_ylabel("Relative abundance"); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout(); return fig


# ---- sample-CSV makers (for the format users upload) ----
def sample_venn_df():
    rng = np.random.default_rng(30)
    allg = np.array([f"GENE{i:04d}" for i in range(220)])
    cols = {}
    for name, n in [("RNAseq", 85), ("Proteomics", 95), ("ATACseq", 75)]:
        members = sorted(rng.choice(allg, n, replace=False))
        cols[name] = members
    m = max(len(v) for v in cols.values())
    return pd.DataFrame({k: v + [None] * (m - len(v)) for k, v in cols.items()})


def sample_enrichment_df():
    rng = np.random.default_rng(32)
    terms = ["Synaptic signaling", "Ion transport", "Immune response", "Cell cycle",
             "Oxidative phosphorylation", "Apoptosis", "MAPK cascade", "Lipid metabolism",
             "DNA repair", "Chemotaxis"]
    return pd.DataFrame({"term": terms,
                         "gene_ratio": np.round(rng.uniform(0.05, 0.45, len(terms)), 3),
                         "count": rng.integers(8, 80, len(terms)),
                         "p_adjust": rng.uniform(1e-6, 3e-2, len(terms))})


def sample_pca_df():
    rng = np.random.default_rng(33)
    groups = np.repeat(["Control", "Disease", "Treated"], 12)
    centres = {"Control": 0, "Disease": 4, "Treated": 2}
    X = np.array([rng.normal(centres[g], 1, 40) + rng.normal(0, 1, 40) for g in groups])
    df = pd.DataFrame(np.round(X, 3), columns=[f"feature{i:02d}" for i in range(40)])
    df.insert(0, "group", groups)
    df.insert(0, "sample", [f"S{i+1:02d}" for i in range(len(groups))])
    return df


def sample_composition_df():
    rng = np.random.default_rng(34)
    taxa = ["Firmicutes", "Bacteroidetes", "Proteobacteria", "Actinobacteria", "Other"]
    comp = rng.uniform(0.2, 1, (10, len(taxa))); comp = comp / comp.sum(1, keepdims=True)
    df = pd.DataFrame(np.round(comp, 4), columns=taxa)
    df.insert(0, "sample", [f"S{i+1}" for i in range(10)])
    return df


FIGURES.update({
    "21 Venn diagram": venn,
    "22 UpSet plot": upset,
    "23 Enrichment dot-plot": enrichment_dotplot,
    "24 PCA / ordination": pca_ordination,
    "25 Stacked composition": stacked_composition,
})


# ============================================================================
# MULTIOMICS real-data versions (upload-driven), inversion-hardened. Each returns
# (figure, warnings). The gallery sample figures are redefined at the end to
# delegate here, so drawing lives in one place.
# ============================================================================

def _draw_upset(sets: dict, preset, width_mm):
    import itertools
    names = list(sets)
    universe = set().union(*sets.values()) if sets else set()
    combos = []
    for r in range(1, len(names) + 1):
        for combo in itertools.combinations(names, r):
            inter = set(universe)
            for n in combo:
                inter &= sets[n]
            for n in names:
                if n not in combo:
                    inter -= sets[n]
            if len(inter):
                combos.append((combo, len(inter)))
    combos.sort(key=lambda x: -x[1]); combos = combos[:14]

    apply_preset(preset) if preset else apply_style()
    fig = plt.figure(figsize=(5.2, 3.4))
    if width_mm:
        set_width_mm(fig, width_mm)
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1.2], hspace=0.05)
    axb = fig.add_subplot(gs[0]); axm = fig.add_subplot(gs[1], sharex=axb)
    xs = np.arange(len(combos))
    axb.bar(xs, [c[1] for c in combos], color=OKABE_ITO[0], width=0.6)
    axb.set_ylabel("Intersection"); axb.spines["bottom"].set_visible(False)
    axb.tick_params(labelbottom=False)
    for j, n in enumerate(names):
        for i, (combo, _) in enumerate(combos):
            axm.scatter(i, j, s=26, color=OKABE_ITO[0] if n in combo else "0.85", zorder=3)
    for i, (combo, _) in enumerate(combos):
        ys = [names.index(n) for n in combo]
        if len(ys) > 1:
            axm.plot([i, i], [min(ys), max(ys)], color=OKABE_ITO[0], lw=1.2, zorder=2)
    axm.set_yticks(range(len(names))); axm.set_yticklabels(names, fontsize=6)
    axm.set_xticks([]); axm.set_ylim(-0.6, len(names) - 0.4)
    for sp in ("top", "right", "bottom"):
        axm.spines[sp].set_visible(False)
    fig.tight_layout(); return fig


def upset_from_data(df, set_cols, preset=None, width_mm=None):
    """UpSet from >=2 columns, each listing the members of one set."""
    warnings = []
    if len(set_cols) < 2:
        raise ValueError("UpSet needs at least 2 set columns.")
    sets = {str(c): set(df[c].dropna().astype(str)) for c in set_cols}
    empty = [k for k, v in sets.items() if not v]
    if empty:
        warnings.append(f"Ignored empty set(s): {', '.join(empty)}.")
        sets = {k: v for k, v in sets.items() if v}
    if len(sets) < 2:
        raise ValueError("Need at least 2 non-empty sets.")
    return _draw_upset(sets, preset, width_mm), warnings


def enrichment_from_data(df, term_col, x_col, size_col=None, color_col=None,
                         preset=None, width_mm=None):
    """GO/KEGG enrichment dot-plot from a results table.

    Inversion guards: x/size/color coerced to numeric; rows missing term or x
    dropped; p-values <=0 floored before -log10; empty -> clear error.
    """
    warnings = []
    d = df[[c for c in [term_col, x_col, size_col, color_col] if c]].copy()
    for c in [x_col, size_col, color_col]:
        if c:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    n0 = len(d); d = d.dropna(subset=[term_col, x_col])
    if len(d) < n0:
        warnings.append(f"Dropped {n0 - len(d)} rows missing term or x-value.")
    if d.empty:
        raise ValueError("No valid enrichment rows after cleaning.")
    d = d.sort_values(x_col).tail(25)
    terms = d[term_col].astype(str).tolist()
    x = d[x_col].to_numpy(float)
    size = (d[size_col].to_numpy(float) if size_col else np.full(len(d), 30.0))
    size = np.nan_to_num(size, nan=np.nanmedian(size) if np.isfinite(np.nanmedian(size)) else 30.0)
    if color_col:
        cv = d[color_col].to_numpy(float)
        if np.nanmin(cv) <= 0:
            warnings.append("Non-positive p.adjust floored before -log10.")
            cv = np.where(cv <= 0, np.nanmin(cv[cv > 0]) if np.any(cv > 0) else 1e-300, cv)
        color = -np.log10(cv); clabel = "-log10 p.adjust"
    else:
        color = x; clabel = x_col
    fig, ax = _new(preset, width_mm, 4.2, max(2.6, 0.28 * len(terms) + 1.2))
    sc = ax.scatter(x, range(len(terms)), s=np.clip(size * 2.5, 10, 300),
                    c=color, cmap="viridis", edgecolor="0.3", linewidth=0.4)
    ax.set_yticks(range(len(terms))); ax.set_yticklabels(terms, fontsize=6)
    ax.set_xlabel(x_col.replace("_", " "))
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02); cb.set_label(clabel, fontsize=6)
    fig.tight_layout(); return fig, warnings


def pca_from_data(df, feature_cols, group_col=None, preset=None, width_mm=None):
    """PCA ordination from a samples x features table.

    Inversion guards: features coerced numeric; constant/all-NaN feature columns
    dropped; rows with any NaN dropped; needs >=2 features and >=3 samples.
    """
    warnings = []
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce")
    const = [c for c in feature_cols if X[c].nunique(dropna=True) <= 1]
    if const:
        warnings.append(f"Dropped {len(const)} constant/empty feature column(s).")
        X = X.drop(columns=const)
    keep = X.dropna()
    if len(keep) < len(X):
        warnings.append(f"Dropped {len(X) - len(keep)} samples with missing values.")
    if X.shape[1] < 2 or len(keep) < 3:
        raise ValueError("PCA needs >=2 numeric features and >=3 complete samples.")
    groups = df.loc[keep.index, group_col].astype(str).to_numpy() if group_col and group_col in df.columns else None
    M = keep.to_numpy(float); M = M - M.mean(0)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    scores = U * S; var = S ** 2 / (S ** 2).sum()
    fig, ax = _new(preset, width_mm, 3.8, 3.2)
    if groups is not None:
        for i, g in enumerate(dict.fromkeys(groups)):
            m = groups == g
            ax.scatter(scores[m, 0], scores[m, 1], s=26, color=OKABE_ITO[i % len(OKABE_ITO)],
                       label=g, edgecolor="white", linewidth=0.5)
        ax.legend(frameon=False, fontsize=6.5)
    else:
        ax.scatter(scores[:, 0], scores[:, 1], s=26, color=OKABE_ITO[0], edgecolor="white", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({var[0]*100:.0f}%)"); ax.set_ylabel(f"PC2 ({var[1]*100:.0f}%)")
    ax.axhline(0, color="0.85", lw=0.6); ax.axvline(0, color="0.85", lw=0.6)
    fig.tight_layout(); return fig, warnings


def stacked_composition_from_data(df, sample_col, category_cols, preset=None, width_mm=None):
    """Stacked relative-abundance bars from a samples x categories table.

    Inversion guards: categories coerced numeric (NaN->0); all-zero samples
    dropped; each sample normalised to sum 1.
    """
    warnings = []
    C = df[category_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    labels = df[sample_col].astype(str).tolist() if sample_col in df.columns else [str(i) for i in range(len(df))]
    tot = C.sum(1)
    zero = (tot <= 0)
    if zero.any():
        warnings.append(f"Dropped {int(zero.sum())} sample(s) with zero total.")
        C = C[~zero.values]; labels = [l for l, z in zip(labels, zero.values) if not z]; tot = tot[~zero.values]
    if C.empty:
        raise ValueError("No samples with non-zero composition.")
    comp = C.to_numpy(float) / tot.to_numpy(float)[:, None]
    fig, ax = _new(preset, width_mm, max(3.4, 0.4 * len(labels) + 1.5), 3.0)
    bottom = np.zeros(len(labels))
    for j, cat in enumerate(category_cols):
        ax.bar(range(len(labels)), comp[:, j], bottom=bottom,
               color=OKABE_ITO[j % len(OKABE_ITO)], label=str(cat), width=0.8)
        bottom += comp[:, j]
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=6, rotation=90 if len(labels) > 12 else 0)
    ax.set_ylabel("Relative abundance"); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.2))
    fig.subplots_adjust(top=0.82)
    return fig, warnings


def sample_upset_df():
    rng = np.random.default_rng(31)
    allg = np.array([f"GENE{i:04d}" for i in range(400)])
    cols = {}
    for name, n in [("RNA", 160), ("Protein", 150), ("ATAC", 140), ("Methyl", 120), ("Metab", 110)]:
        cols[name] = sorted(rng.choice(allg, n, replace=False))
    m = max(len(v) for v in cols.values())
    return pd.DataFrame({k: v + [None] * (m - len(v)) for k, v in cols.items()})


# ---- redefine gallery sample figures to delegate to the from_data drawing ----
def venn(preset=None, width_mm=None):
    return venn_from_data(sample_venn_df(), ["RNAseq", "Proteomics", "ATACseq"],
                          preset=preset, width_mm=width_mm)[0]


def upset(preset=None, width_mm=None):
    return upset_from_data(sample_upset_df(), ["RNA", "Protein", "ATAC", "Methyl", "Metab"],
                           preset=preset, width_mm=width_mm)[0]


def enrichment_dotplot(preset=None, width_mm=None):
    return enrichment_from_data(sample_enrichment_df(), "term", "gene_ratio", "count",
                                "p_adjust", preset=preset, width_mm=width_mm)[0]


def pca_ordination(preset=None, width_mm=None):
    df = sample_pca_df()
    feats = [c for c in df.columns if c.startswith("feature")]
    return pca_from_data(df, feats, "group", preset=preset, width_mm=width_mm)[0]


def stacked_composition(preset=None, width_mm=None):
    df = sample_composition_df()
    cats = [c for c in df.columns if c != "sample"]
    return stacked_composition_from_data(df, "sample", cats, preset=preset, width_mm=width_mm)[0]


FIGURES.update({
    "21 Venn diagram": venn, "22 UpSet plot": upset, "23 Enrichment dot-plot": enrichment_dotplot,
    "24 PCA / ordination": pca_ordination, "25 Stacked composition": stacked_composition,
})
REAL_DATA_FIGURES.update({
    "21 Venn diagram": "venn", "22 UpSet plot": "upset", "23 Enrichment dot-plot": "enrichment",
    "24 PCA / ordination": "pca", "25 Stacked composition": "composition",
})


# ============================================================================
# CLINICAL / BIOMEDICAL figures: ridgeline, chord, Kaplan-Meier (+log-rank), forest
# ============================================================================

def ridgeline(preset=None, width_mm=None):
    """Ridgeline (joyplot): stacked density distributions across an ordered factor."""
    fig, ax = _new(preset, width_mm, 3.8, 3.6)
    rng = _rng(40)
    groups = ["Day 1", "Day 3", "Day 7", "Day 14", "Day 21", "Day 28"]
    for i, g in enumerate(groups):
        d = rng.normal(i * 0.55, 1.0 + 0.08 * i, 220)
        kde = stats.gaussian_kde(d)
        xs = np.linspace(d.min(), d.max(), 200); ys = kde(xs); ys = ys / ys.max() * 1.5
        base = len(groups) - i
        ax.fill_between(xs, base, base + ys, color=OKABE_ITO[i % len(OKABE_ITO)],
                        alpha=0.75, lw=0.8, edgecolor="white", zorder=len(groups) - i)
    ax.set_yticks([len(groups) - i for i in range(len(groups))])
    ax.set_yticklabels(groups, fontsize=6)
    ax.set_xlabel("Expression"); ax.set_ylabel("")
    for sp in ("left", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); return fig


def chord_diagram(preset=None, width_mm=None):
    """Chord diagram: weighted flows between categories, arranged on a circle."""
    fig, ax = _new(preset, width_mm, 3.6, 3.6)
    rng = _rng(41)
    n = 6; labels = [f"Cell{i+1}" for i in range(n)]
    flow = rng.uniform(0, 1, (n, n)); flow = (flow + flow.T) / 2; np.fill_diagonal(flow, 0)
    totals = flow.sum(1); grand = totals.sum(); gap = 0.03 * 2 * np.pi
    spans = totals / grand * (2 * np.pi - n * gap)
    start = {}; a = np.pi / 2
    for i in range(n):
        start[i] = a; a -= spans[i] + gap
    for i in range(n):
        ax.add_patch(Wedge((0, 0), 1.0, np.degrees(start[i] - spans[i]), np.degrees(start[i]),
                           width=0.08, facecolor=OKABE_ITO[i % len(OKABE_ITO)], edgecolor="white"))
        mid = start[i] - spans[i] / 2
        ax.text(1.12 * np.cos(mid), 1.12 * np.sin(mid), labels[i], ha="center", va="center", fontsize=6)
    pt = lambda ang: (0.9 * np.cos(ang), 0.9 * np.sin(ang))
    sub = {i: start[i] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            f = flow[i, j]
            if f / grand < 0.01:
                continue
            w = f / grand * (2 * np.pi - n * gap)
            ai0 = sub[i]; ai1 = sub[i] - w; sub[i] = ai1
            aj0 = sub[j]; aj1 = sub[j] - w; sub[j] = aj1
            verts = [pt(ai0), (0, 0), pt(aj0), pt(aj1), (0, 0), pt(ai1), pt(ai0)]
            codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO,
                     Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]
            ax.add_patch(mpatches.PathPatch(Path(verts, codes),
                         facecolor=OKABE_ITO[i % len(OKABE_ITO)], edgecolor="none", alpha=0.4))
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout(); return fig


# ---- Kaplan-Meier survival + log-rank (real analysis) ----
def _km_estimate(dur, ev):
    dur = np.asarray(dur, float); ev = np.asarray(ev, int)
    times = np.array(sorted(set(dur)))
    t_out = [0.0]; s_out = [1.0]; s = 1.0; cens = []
    for t in times:
        n = int((dur >= t).sum())
        d = int(((dur == t) & (ev == 1)).sum())
        c = int(((dur == t) & (ev == 0)).sum())
        if n > 0 and d > 0:
            s *= (1 - d / n)
        t_out.append(float(t)); s_out.append(s)
        if c > 0:
            cens.append((float(t), s))
    return np.array(t_out), np.array(s_out), cens


def _logrank(dur, ev, grp):
    dur = np.asarray(dur, float); ev = np.asarray(ev, int); grp = np.asarray(grp)
    g1 = grp == np.unique(grp)[0]
    O1 = E1 = V = 0.0
    for t in np.array(sorted(set(dur[ev == 1]))):
        at = dur >= t
        n = int(at.sum()); n1 = int((at & g1).sum())
        d = int(((dur == t) & (ev == 1)).sum()); d1 = int(((dur == t) & (ev == 1) & g1).sum())
        if n > 1:
            E1 += d * n1 / n
            V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
        O1 += d1
    chi2 = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return chi2, float(stats.chi2.sf(chi2, 1))


def km_from_data(df, time_col, event_col, group_col=None, preset=None, width_mm=None):
    """Kaplan-Meier survival curves (+ log-rank p for two groups).

    Guards: time/event coerced numeric; event must be 0/1 (censored/event);
    rows with missing time or event dropped.
    """
    warnings = []
    d = df[[c for c in [time_col, event_col, group_col] if c]].copy()
    d[time_col] = pd.to_numeric(d[time_col], errors="coerce")
    d[event_col] = pd.to_numeric(d[event_col], errors="coerce")
    n0 = len(d); d = d.dropna(subset=[time_col, event_col])
    if len(d) < n0:
        warnings.append(f"Dropped {n0 - len(d)} rows missing time/event.")
    ev = d[event_col].to_numpy()
    if not set(np.unique(ev)).issubset({0, 1}):
        warnings.append("Event column should be 0 (censored) / 1 (event); "
                        "non-zero values treated as events.")
        d[event_col] = (d[event_col] != 0).astype(int)
    fig, ax = _new(preset, width_mm, 3.8, 3.2)
    groups = list(dict.fromkeys(d[group_col].astype(str))) if group_col and group_col in d.columns else [None]
    for i, g in enumerate(groups):
        sub = d if g is None else d[d[group_col].astype(str) == g]
        t, s, cens = _km_estimate(sub[time_col].to_numpy(), sub[event_col].to_numpy())
        c = OKABE_ITO[i % len(OKABE_ITO)]
        ax.step(t, s, where="post", color=c, lw=1.4, label=(g or "All"))
        if cens:
            ct, cs = zip(*cens)
            ax.scatter(ct, cs, marker="|", s=30, color=c, zorder=3)
    ax.set_xlabel(time_col.replace("_", " ")); ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    cap = ""
    if group_col and len(groups) == 2:
        chi2, p = _logrank(d[time_col].to_numpy(), d[event_col].to_numpy(),
                           d[group_col].astype(str).to_numpy())
        from neurofig_core import fmt_p
        cap = f"log-rank: χ²={chi2:.2f}, p {fmt_p(p)}"
    if groups != [None]:
        ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    if cap:
        fig.text(0.5, -0.04, cap, ha="center", va="top", fontsize=6.4, color="0.25")
    fig.tight_layout(); return fig, warnings


def kaplan_meier(preset=None, width_mm=None):
    return km_from_data(sample_km_df(), "time_months", "event", "arm",
                        preset=preset, width_mm=width_mm)[0]


# ---- Forest plot ----
def forest_from_data(df, label_col, estimate_col, low_col, high_col,
                     ref=0.0, xlabel="Effect (95% CI)", preset=None, width_mm=None):
    """Forest plot of effect estimates with 95% CIs (meta-analysis / regression)."""
    warnings = []
    d = df[[label_col, estimate_col, low_col, high_col]].copy()
    for c in [estimate_col, low_col, high_col]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    n0 = len(d); d = d.dropna(subset=[estimate_col, low_col, high_col])
    if len(d) < n0:
        warnings.append(f"Dropped {n0 - len(d)} rows with missing estimate/CI.")
    if d.empty:
        raise ValueError("No valid rows for the forest plot.")
    labels = d[label_col].astype(str).tolist()
    est = d[estimate_col].to_numpy(float); lo = d[low_col].to_numpy(float); hi = d[high_col].to_numpy(float)
    y = np.arange(len(labels))[::-1]
    fig, ax = _new(preset, width_mm, 4.0, max(2.4, 0.32 * len(labels) + 1.0))
    ax.errorbar(est, y, xerr=[est - lo, hi - est], fmt="s", ms=5, color=OKABE_ITO[0],
                ecolor="0.4", lw=1.0, capsize=2)
    ax.axvline(ref, color="0.5", ls="--", lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel(xlabel); ax.set_ylim(-0.6, len(labels) - 0.4)
    fig.tight_layout(); return fig, warnings


def forest(preset=None, width_mm=None):
    return forest_from_data(sample_forest_df(), "study", "estimate", "ci_low", "ci_high",
                            ref=0.0, preset=preset, width_mm=width_mm)[0]


# ---- sample makers ----
def sample_km_df():
    rng = np.random.default_rng(42)
    rows = []
    for arm, scale, cens in [("Standard", 12, 0.30), ("Experimental", 22, 0.35)]:
        for _ in range(60):
            tt = rng.exponential(scale)
            fu = rng.uniform(6, 40)  # follow-up / censoring time
            time = min(tt, fu); event = int(tt <= fu and rng.random() > cens)
            rows.append({"arm": arm, "time_months": round(time, 1), "event": event})
    return pd.DataFrame(rows)


def sample_forest_df():
    rng = np.random.default_rng(43)
    studies = ["Anderson 2019", "Brown 2020", "Chen 2020", "Dubois 2021",
               "Evans 2021", "Farah 2022", "Gupta 2023"]
    est = rng.normal(0.35, 0.25, len(studies))
    se = rng.uniform(0.1, 0.3, len(studies))
    return pd.DataFrame({"study": studies, "estimate": np.round(est, 3),
                         "ci_low": np.round(est - 1.96 * se, 3),
                         "ci_high": np.round(est + 1.96 * se, 3)})


FIGURES.update({
    "26 Ridgeline": ridgeline, "27 Chord diagram": chord_diagram,
    "28 Kaplan-Meier survival": kaplan_meier, "29 Forest plot": forest,
})
REAL_DATA_FIGURES.update({"28 Kaplan-Meier survival": "km", "29 Forest plot": "forest"})
