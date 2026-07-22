"""
standardcurve.py — fit a standard curve and back-calculate unknown samples.

Run known-concentration standards through an assay, measure the signal, fit a
curve, then read unknown concentrations off it. The workhorse of every wet lab
(protein assays, ELISA, qPCR), usually done in fragile Excel sheets.

Three models cover the common assays:
  • "linear"  — signal ∝ concentration (BCA / Bradford protein assays)
  • "qpcr"    — Ct vs log10(quantity) is linear; also reports amplification
                efficiency E = 10^(-1/slope) − 1
  • "4pl"     — sigmoidal (ELISA), the four-parameter logistic

The point of a standard curve is the *inverse*: given an unknown's signal, return
its concentration. Every model implements back_calculate() and is validated by a
round-trip against known values (see tests/test_standardcurve.py).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from neurofig_core import OKABE_ITO, apply_preset, apply_style, set_width_mm


def _fourpl(logx, bottom, top, logic50, hill):
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((logic50 - logx) * hill))


@dataclass
class StandardCurveFit:
    model: str
    params: dict
    r_squared: float
    efficiency: float | None = None          # qPCR only
    _popt: tuple = field(default=(), repr=False)

    def predict(self, x):
        """Signal expected at concentration(s) x."""
        x = np.asarray(x, float)
        if self.model == "linear":
            return self.params["slope"] * x + self.params["intercept"]
        if self.model == "qpcr":
            return self.params["slope"] * np.log10(x) + self.params["intercept"]
        return _fourpl(np.log10(x), *self._popt)

    def back_calculate(self, y):
        """Concentration(s) for measured signal(s) y (the inverse of the curve).

        Returns NaN where y falls outside the model's invertible range (e.g. a
        4PL signal above the top plateau), rather than a fabricated number.
        """
        y = np.asarray(y, float)
        if self.model == "linear":
            m, b = self.params["slope"], self.params["intercept"]
            return (y - b) / m
        if self.model == "qpcr":
            m, b = self.params["slope"], self.params["intercept"]
            return 10.0 ** ((y - b) / m)
        bottom, top, logic50, hill = self._popt
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (top - bottom) / (y - bottom) - 1.0
            logx = logic50 - np.log10(ratio) / hill
            x = 10.0 ** logx
            x = np.where(ratio > 0, x, np.nan)   # y must lie between plateaus
        return x


def _r_squared(y, yhat):
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def fit_standard_curve(conc, signal, model: str = "linear") -> StandardCurveFit:
    """Fit standards. model = 'linear' | 'qpcr' | '4pl'.

    conc, signal : equal-length arrays of known concentrations and measured
                   signals (replicates allowed). qpcr/4pl require conc > 0.
    """
    conc = np.asarray(conc, float)
    signal = np.asarray(signal, float)
    m = np.isfinite(conc) & np.isfinite(signal)
    if model in ("qpcr", "4pl"):
        m &= conc > 0
    conc, signal = conc[m], signal[m]
    if conc.size < (4 if model == "4pl" else 2):
        raise ValueError(f"Not enough valid points for a {model} fit.")

    if model == "linear":
        slope, intercept = np.polyfit(conc, signal, 1)
        fit = StandardCurveFit("linear", {"slope": float(slope), "intercept": float(intercept)},
                               0.0)
        fit.r_squared = _r_squared(signal, fit.predict(conc))
        return fit

    if model == "qpcr":
        logx = np.log10(conc)
        slope, intercept = np.polyfit(logx, signal, 1)
        eff = 10.0 ** (-1.0 / slope) - 1.0 if slope != 0 else float("nan")
        fit = StandardCurveFit("qpcr", {"slope": float(slope), "intercept": float(intercept)},
                               0.0, efficiency=float(eff))
        fit.r_squared = _r_squared(signal, fit.predict(conc))
        return fit

    if model == "4pl":
        logx = np.log10(conc)
        b0, t0 = float(np.min(signal)), float(np.max(signal))
        sign = 1.0 if np.corrcoef(logx, signal)[0, 1] >= 0 else -1.0
        popt = None
        for guess in ([b0, t0, float(np.median(logx)), sign],
                      [b0, t0, float(np.median(logx)), -sign]):
            try:
                popt, _ = curve_fit(_fourpl, logx, signal, p0=guess, maxfev=20000)
                break
            except (RuntimeError, ValueError):
                continue
        if popt is None:
            raise RuntimeError("4PL standard-curve fit did not converge.")
        fit = StandardCurveFit("4pl",
                               {"bottom": float(popt[0]), "top": float(popt[1]),
                                "ic50": float(10 ** popt[2]), "hill": float(popt[3])},
                               0.0, _popt=tuple(popt))
        fit.r_squared = _r_squared(signal, _fourpl(logx, *popt))
        return fit

    raise ValueError("model must be 'linear', 'qpcr', or '4pl'.")


def make_standard_curve_figure(conc, signal, fit: StandardCurveFit,
                               unknowns: list | None = None, *,
                               xlabel: str = "Concentration", ylabel: str = "Signal",
                               title: str = "", preset: str | None = None,
                               width_mm: float | None = None) -> plt.Figure:
    """Standards + fitted curve; optionally mark back-calculated unknown samples."""
    apply_preset(preset) if preset else apply_style()
    conc = np.asarray(conc, float)
    signal = np.asarray(signal, float)
    logscale = fit.model in ("qpcr", "4pl")
    c0, c1 = OKABE_ITO[0], OKABE_ITO[1]

    fig, ax = plt.subplots(figsize=(3.4, 2.8))
    if width_mm:
        set_width_mm(fig, width_mm)
    if logscale:
        ax.set_xscale("log")

    ax.scatter(conc, signal, s=18, facecolor="white", edgecolor=c0, linewidth=1.0,
               zorder=3, label="Standards")
    lo, hi = np.min(conc[conc > 0]) if logscale else np.min(conc), np.max(conc)
    xx = np.logspace(np.log10(lo), np.log10(hi), 200) if logscale else np.linspace(lo, hi, 200)
    ax.plot(xx, fit.predict(xx), color=c0, lw=1.3, zorder=2)

    if unknowns:
        ys = np.array([u[1] if isinstance(u, (tuple, list)) else u for u in unknowns], float)
        xs = fit.back_calculate(ys)
        ax.scatter(xs, ys, s=26, marker="D", color=c1, zorder=4, label="Unknowns")

    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=6.5, loc="best")

    eq = {"linear": f"y = {fit.params.get('slope', 0):.3g}·x + {fit.params.get('intercept', 0):.3g}",
          "qpcr": f"Ct = {fit.params.get('slope', 0):.3g}·log10(x) + {fit.params.get('intercept', 0):.3g}",
          "4pl": "4-parameter logistic"}[fit.model]
    cap = f"{eq}   ·   R² = {fit.r_squared:.4f}"
    if fit.efficiency is not None:
        cap += f"   ·   efficiency = {fit.efficiency*100:.0f}%"
    fig.text(0.5, -0.04, cap, ha="center", va="top", fontsize=6.2, color="0.25")
    fig.tight_layout()
    return fig
