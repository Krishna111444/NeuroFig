# NeuroFig — Validation & Reproducibility Dossier

*Scientific validation of every computational component. Last validated 2026-07-23,
Python 3.13.13, all 66 tests passing.*

This document exists so a reviewer, collaborator, or paying user can verify — not
take on faith — that NeuroFig computes the right answer, does so reproducibly, and
refuses to produce a wrong one. It doubles as the methods/validation reference for
any manuscript that uses the tool.

---

## 1. Validation philosophy (first principles)

A scientific tool earns trust on four properties. We hold each one and show the
evidence:

| Property | What it means | How we guarantee it |
|---|---|---|
| **Correctness** | Methods implement established science correctly | Every method validated against **ground truth** (known parameters recovered, or agreement with a reference implementation) |
| **Reproducibility** | Same input → same output, forever | Deterministic code, fixed seeds, recorded library versions, open source, no network/AI in the compute path |
| **Transparency** | The user always knows what was done | The test used is named in every result; assumptions and limitations are stated, not hidden |
| **Robustness** | Messy/adversarial input degrades gracefully | Failure modes enumerated by **inversion** and each one guarded + tested |

**Architecture that enforces this:** all statistics and drawing live in a pure,
deterministic core (`neurofig_core.py`, `doseresponse.py`, `standardcurve.py`,
`cluster_stats.py`, `gasteiger.py`, `labcalc.py`, `pdbqt_tools.py`). There is **no
machine-learning / no LLM anywhere in the compute path** — nothing is inferred or
"generated," so nothing can be hallucinated. The UI (`app.py`) only calls the core.

---

## 2. Component-by-component validation

Each row states the claim, the **ground truth** it was checked against, the result,
and the test that re-verifies it.

### Statistics — group comparison (`neurofig_core.py`)
- **Test auto-selection:** 2 groups → Welch's *t* (never assumes equal variance);
  3+ → one-way ANOVA + Tukey HSD; non-normal (Shapiro) → Mann-Whitney / Kruskal-
  Wallis + Dunn. The test actually used is always reported to the user.
- **Dunn's post-hoc (own implementation):** validated **atom-for-atom against
  `scikit-posthocs` to < 1e-6** on Holm-adjusted p-values.
- **p-value reporting:** `fmt_p` never prints `p = 0` (invalid); values < 0.001 →
  "< 0.001". Verified `fmt_p(0.0)` and `fmt_p(1e-300)` both give "< 0.001".
- Tests: `test_core.py` (test_infer_and_three_group_anova, test_nonparametric_*,
  test_dunn_matches_kruskal_direction, test_fmt_p_never_zero).

### Dose-response IC50/EC50 — 4-parameter logistic (`doseresponse.py`)
- Same variable-slope 4PL model as GraphPad Prism.
- **Ground truth — parameter recovery:** on noise-free simulated data the fit
  recovers the true IC50 and Hill slope **to 4 significant figures (R² = 1.0000)**;
  on noisy data, IC50 within ~5% of truth. Agonist curves auto-label **EC50**,
  inhibitors **IC50**; the 95% CI brackets the estimate.
- Tests: `test_doseresponse.py` (recovers_agonist_ec50, recovers_inhibitor_ic50,
  noise_free_is_exact, ci_and_predict_and_guards).

### Standard curve + back-calculation (`standardcurve.py`)
- Models: linear (protein assays), qPCR (Ct vs log-quantity), 4PL (ELISA).
- **Ground truth:** linear slope recovered exactly; **qPCR efficiency = 100.0%**
  from a slope of −3.323 (true −3.32); 4PL back-calculation round-trips a known
  concentration to **5.0000**. A signal outside the curve's range returns **NaN,
  never a fabricated concentration.**
- Tests: `test_standardcurve.py` (linear_recovery_and_roundtrip, qpcr_efficiency,
  4pl_roundtrip_and_out_of_range).

### Cluster-based permutation test (`cluster_stats.py`)
- Maris–Oostenveld (2007) method for comparing two time-courses with family-wise
  error control across time.
- **Ground truth — the defining property of a valid test:** under the null (two
  conditions from the *same* distribution), the false-positive rate is **~0.06**
  against a nominal 0.05. With a real effect injected, it recovers and correctly
  **localizes** the significant window (p < 0.01).
- Tests: `test_cluster.py` (cluster_finder_exact, power_detects_and_localizes,
  false_positive_rate_controlled).

### Gasteiger partial charges (`gasteiger.py`)
- From-scratch PEOE implementation (Gasteiger–Marsili 1980) **+ carbonyl-
  conjugation (resonance) handling.**
- **Ground truth — Open Babel 3.1.0:** exact agreement (**0.0000 e**) on 12 of 13
  test molecules including ibuprofen, aspirin (two carbonyls), glycine, esters.
  Primary amide residual **0.04 e** (within Gasteiger's own approximation, and
  documented). Non-conjugated molecules unchanged — no false positives.
- Tests: `test_gasteiger.py` (methane/water regression, convergence,
  matches_openbabel_on_nonresonant, resonance_matches_openbabel_on_carboxyls).

### Bench calculators (`labcalc.py`)
- 13 deterministic calculations across solution prep, nucleic acids, protein,
  spectrophotometry, buffers, microbiology, pharmacology.
- **Ground truth — known reference values:** glycine MW = 75.07 Da, diglycine =
  132.12 Da, qPCR copy number = 9.26×10⁸, FDA human-equivalent-dose Km table,
  Beer–Lambert, Henderson–Hasselbalch, doubling time. All match.
- Tests: `test_labcalc.py` (solution, nucleic_acids, protein, spectro_and_buffer,
  micro_and_pharma).

### Structure prep — PDB → PDBQT (`pdbqt_tools.py`)
- Pure-Python PDB cleaning (waters, chains, alternate locations, hetero, ligand
  extraction) + Open Babel conversion.
- **Ground truth — real conversion:** ibuprofen converts to a valid PDBQT with a
  torsion tree, AutoDock atom types, **and non-zero Gasteiger charges** (a real
  bug — all-zero charges that silently break docking — was found and fixed here
  once the backend was live).
- Tests: `test_pdbqt.py` (summary, clean_*, extract_ligand, guard when missing,
  real_conversion_assigns_charges).

### Figures (`neurofig_core.py`, `journal_figures.py`)
- **Colour:** Okabe-Ito colourblind-safe palette throughout.
- **Honesty of representation:** bars are mean±SEM with individual points shown;
  significance brackets are placed from the **actual data extent** (never overlap
  data); signed data (ΔF/F, z-scores) get a correct baseline.
- **Exact journal sizing — validated by measurement:** a Nature single-column
  export measures **89.00 mm** (PDF MediaBox), not "about right."
- **Real-data figures (volcano, raincloud, correlation heatmap)** are inversion-
  hardened (Section 3). The **20-figure gallery** is rendering-tested and export-
  tested at journal presets.
- Tests: `test_core.py` (figure exports, journal width exact), `test_journal_figures.py`,
  `test_realdata.py`.

---

## 3. Inversion audit — "how could this be wrong?"

For every user-facing computation we enumerated the ways it could crash or silently
mislead, then guarded and **tested** each. Selected cases:

| Component | Failure mode | Guard (all tested) |
|---|---|---|
| Volcano | `p = 0` → −log10 = ∞ | floored to smallest non-zero p, with a warning |
| Volcano | p outside [0,1]; linear FC ≤ 0 | dropped as invalid |
| Volcano | user's p already −log10 | option prevents a double transform |
| Raincloud | group n<2 or zero variance | falls back to points (no KDE singular crash) |
| Correlation | constant column (r undefined) | dropped with a warning, not silent NaN |
| Correlation | missing values | pandas pairwise-complete correlation |
| Group stats | single group / empty group | clear error / dropped with warning |
| Paired stats | < 2 complete pairs | clear error |
| Standard curve | signal outside curve range | NaN, never an extrapolated fake |
| p-values | extreme significance clamped to 0 | reported "< 0.001", never "p = 0" |
| Time-course | AUC window < 2 points | clear error |

Inversion tests: `test_realdata.py` (14 tests feeding the pathological inputs) plus
guard tests across the other suites.

---

## 4. Reproducibility

- **Determinism (verified):** repeated runs of the group comparison, dose-response
  fit, Gasteiger charges, seeded cluster test, calculators, and figure statistics
  each produce **identical** output. Same input → same answer, always.
- **Seeds:** all sample-data generators and the permutation test take explicit
  seeds; nothing depends on wall-clock time.
- **Recorded environment:** exact library versions in `requirements-lock.txt`
  (numpy 2.4.6, scipy 1.18.0, pandas 2.3.3, matplotlib 3.11.0, statsmodels 0.14.6,
  streamlit 1.60.0; Python 3.13.13). Open Babel 3.1.0 for structure prep.
- **Open source & offline:** the entire compute path is inspectable Python and runs
  with no network call and no external model — a reviewer can re-run everything.

Reproduce this whole dossier:
```bash
pip install -r requirements-lock.txt
python -m tests.test_core && python -m tests.test_doseresponse && \
python -m tests.test_standardcurve && python -m tests.test_cluster && \
python -m tests.test_gasteiger && python -m tests.test_labcalc && \
python -m tests.test_pdbqt && python -m tests.test_journal_figures && \
python -m tests.test_realdata
```

---

## 5. Honest limitations & scope

Stated plainly, because hiding them is what makes a tool *un*scientific:

1. **Unequal variance, 3+ groups.** We flag heteroscedasticity with a Levene test
   but still run ANOVA + Tukey. The strictly-correct choice (Welch's ANOVA +
   Games-Howell) is **not yet implemented** — it is flagged, not silently ignored.
2. **Normality decision.** Shapiro–Wilk drives the parametric/non-parametric default
   only; it is under-powered at small *n* and over-sensitive at large *n*, so the
   user can always override the choice.
3. **Primary amides** in Gasteiger charges differ from Open Babel by ~0.04 e (within
   the method's own approximation).
4. **Volcano multiplicity.** We plot the p-values the user supplies; we do **not**
   compute FDR/adjusted p — supply adjusted p (e.g. `padj`) if you need it.
5. **The 20-figure gallery is visualization, not inference.** Seventeen of the
   twenty are rendering tools shown with sample data; only the three real-data
   figures (volcano, raincloud, correlation heatmap) validate and analyse uploads.
   The gallery figures must not be read as statistical claims.
6. **No AI layer.** The planned plain-English interface is not built — which is why
   there is currently *zero* hallucination surface. If it is added, it will only
   translate English → a structured request that the validated core executes.

---

## 6. Working status (2026-07-23)

| Area | Status |
|---|---|
| Deterministic engine (stats, fits, charges, calculators) | ✅ validated, 66 tests |
| 4 analysis workflows (group, time-course+cluster, IC50/EC50, standard curve) | ✅ in-app, validated |
| 20-figure journal gallery (+3 real-data) | ✅ in-app, render/export tested |
| PDB→PDBQT (clean + convert) | ✅ validated with live Open Babel |
| Journal presets + exact-mm vector export | ✅ measured exact |
| Paywall (watermark + offline signed licenses) | ✅ built & tested |
| Self-hosting runbook | ✅ `SELF_HOSTING.md` |
| App runs end-to-end in a browser | ✅ verified live |

---

## 7. Test inventory
`test_core` 23 · `test_pdbqt` 6 · `test_cluster` 3 · `test_labcalc` 5 ·
`test_doseresponse` 4 · `test_standardcurve` 3 · `test_gasteiger` 5 ·
`test_journal_figures` 3 · `test_realdata` 14 — **66 total, 0 failures.**
