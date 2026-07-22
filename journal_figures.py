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
