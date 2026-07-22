"""
doseresponse.py — IC50 / EC50 dose-response analysis (4-parameter logistic).

Fits the industry-standard variable-slope four-parameter logistic (the same model
as GraphPad Prism's "log(agonist/inhibitor) vs response — variable slope"):

    Y = Bottom + (Top - Bottom) / (1 + 10^((logIC50 - X)·HillSlope))     X = log10[conc]

and returns IC50/EC50 with a 95% confidence interval, the Hill slope, the
plateaus, and R². Direction (agonist EC50 vs inhibitor IC50) falls out of the fit.

This is a real algorithm, validated by recovering known parameters from simulated
data (see tests/test_doseresponse.py) — not a wrapper and not a guess.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from neurofig_core import OKABE_ITO, apply_preset, apply_style, set_width_mm


def _fourpl(logx, bottom, top, logic50, hill):
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((logic50 - logx) * hill))


@dataclass
class DoseResponseFit:
    ic50: float                       # = EC50 for an agonist; the inflection conc
    ic50_ci: tuple[float, float]      # 95% CI (asymmetric in linear space)
    hill: float                       # Hill slope (signed)
    top: float
    bottom: float
    r_squared: float
    direction: str                    # "increasing" (agonist) | "decreasing" (inhibitor)
    label: str                        # "EC50" or "IC50"
    n: int
    popt: tuple

    def predict(self, conc):
        """Fitted response at concentration(s) `conc` (linear units, > 0)."""
        return _fourpl(np.log10(np.asarray(conc, float)), *self.popt)


def fit_dose_response(conc, response, kind: str = "auto") -> DoseResponseFit:
    """Fit the 4-parameter logistic to (concentration, response) data.

    conc, response : equal-length arrays; multiple rows at the same concentration
                     (replicates) are allowed and all used in the fit.
    kind           : 'auto' labels IC50 vs EC50 from the fitted direction;
                     'ic50'/'inhibitor' or 'ec50'/'agonist' force the label.
    """
    conc = np.asarray(conc, float)
    response = np.asarray(response, float)
    m = np.isfinite(conc) & np.isfinite(response) & (conc > 0)
    conc, response = conc[m], response[m]
    if conc.size < 4:
        raise ValueError("Need at least 4 valid (concentration>0, response) points.")

    logx = np.log10(conc)
    b0, t0 = float(np.min(response)), float(np.max(response))
    slope_sign = 1.0 if np.corrcoef(logx, response)[0, 1] >= 0 else -1.0
    p0 = [b0, t0, float(np.median(logx)), slope_sign]

    popt = None
    for guess in (p0, [b0, t0, p0[2], -slope_sign]):
        try:
            popt, pcov = curve_fit(_fourpl, logx, response, p0=guess, maxfev=20000)
            break
        except (RuntimeError, ValueError):
            continue
    if popt is None:
        raise RuntimeError("Dose-response fit did not converge; check the data.")

    bottom, top, logic50, hill = (float(v) for v in popt)
    se_log = float(np.sqrt(pcov[2, 2])) if np.all(np.isfinite(pcov)) else float("nan")
    ic50 = 10.0 ** logic50
    if np.isfinite(se_log):
        ci = (10.0 ** (logic50 - 1.96 * se_log), 10.0 ** (logic50 + 1.96 * se_log))
    else:
        ci = (float("nan"), float("nan"))

    yhat = _fourpl(logx, *popt)
    ss_res = float(np.sum((response - yhat) ** 2))
    ss_tot = float(np.sum((response - np.mean(response)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    lo, hi = _fourpl(logx.min() - 3, *popt), _fourpl(logx.max() + 3, *popt)
    increasing = hi >= lo
    direction = "increasing" if increasing else "decreasing"
    if kind == "auto":
        label = "EC50" if increasing else "IC50"
    elif kind in ("ec50", "agonist"):
        label = "EC50"
    elif kind in ("ic50", "inhibitor"):
        label = "IC50"
    else:
        raise ValueError("kind must be auto/ec50/agonist/ic50/inhibitor")

    return DoseResponseFit(ic50=ic50, ic50_ci=ci, hill=hill, top=top, bottom=bottom,
                           r_squared=r2, direction=direction, label=label,
                           n=int(conc.size), popt=tuple(popt))


def make_dose_response_figure(conc, response, fit: DoseResponseFit, *,
                              ylabel: str = "Response (%)",
                              xlabel: str = "Concentration",
                              title: str = "",
                              preset: str | None = None,
                              width_mm: float | None = None) -> plt.Figure:
    """Semi-log dose-response: data points, fitted sigmoid, and the IC50/EC50 mark."""
    apply_preset(preset) if preset else apply_style()
    conc = np.asarray(conc, float)
    response = np.asarray(response, float)
    m = np.isfinite(conc) & np.isfinite(response) & (conc > 0)
    conc, response = conc[m], response[m]
    c = OKABE_ITO[0]

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    if width_mm:
        set_width_mm(fig, width_mm)
    ax.set_xscale("log")

    # mean ± SEM per unique concentration (points), plus the raw fit
    uniq = np.unique(conc)
    means = np.array([response[conc == u].mean() for u in uniq])
    sems = np.array([response[conc == u].std(ddof=1) / np.sqrt(np.sum(conc == u))
                     if np.sum(conc == u) > 1 else 0.0 for u in uniq])
    ax.errorbar(uniq, means, yerr=sems, fmt="o", ms=5, color=c, capsize=3,
                lw=1.0, mfc="white", mec=c, zorder=3)

    xx = np.logspace(np.log10(conc.min()), np.log10(conc.max()), 200)
    ax.plot(xx, fit.predict(xx), color=c, lw=1.4, zorder=2)

    ax.axvline(fit.ic50, color="0.5", ls="--", lw=0.8, zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")

    ci_txt = ""
    if np.isfinite(fit.ic50_ci[0]):
        ci_txt = f"  (95% CI {fit.ic50_ci[0]:.3g}–{fit.ic50_ci[1]:.3g})"
    fig.text(0.5, -0.04,
             f"{fit.label} = {fit.ic50:.3g}{ci_txt}   ·   Hill = {fit.hill:.2f}   ·   R² = {fit.r_squared:.3f}",
             ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig
