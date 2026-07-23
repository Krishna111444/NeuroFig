# NeuroFig

**Upload your data → get a journal-ready figure with the correct statistics.**
Colourblind-safe, exact-mm vector export, no coding. A validated, deterministic
figure-and-analysis tool for the life sciences.

The statistics and drawing are done by **deterministic, tested code** — there is no
machine-learning in the compute path, so nothing can be hallucinated. Every method
is validated against ground truth. See **[VALIDATION.md](VALIDATION.md)**.

---

## What it does

**Analysis workflows** (each computes the correct test and draws the figure):
- **Group comparison** — auto-selects Welch's *t* / ANOVA+Tukey / Mann-Whitney /
  Kruskal-Wallis+Dunn from a normality check, and tells you which it used.
- **Time-course (ΔF/F)** — mean±SEM traces + a **cluster-based permutation test**
  (validated false-positive control).
- **Dose-response (IC50/EC50)** — 4-parameter logistic, Prism-grade, with CI.
- **Standard curve** — linear / qPCR (with efficiency) / 4PL, with back-calculation
  of unknowns.

**Figure gallery** — the 20 recurring top-journal figure types (volcano, Manhattan,
raincloud, swimmer, Circos, Sankey, treemap, split violin, …). Volcano, raincloud,
and correlation heatmap accept your own uploads.

**Structure prep** — clean a PDB and convert to **PDBQT** for AutoDock/Vina (via
Open Babel), with correct Gasteiger charges.

**Bench calculators** — 13 deterministic, reference-checked calculations (molarity,
dilutions, DNA copy number, primer/protein properties, Beer–Lambert, CFU, doubling
time, FDA human-equivalent dose, …).

**Journal presets** (Nature/Cell/eLife/Science) with exact figure-width-in-mm export,
in-figure editing, and a watermark/licence paywall.

---

## Run it locally
```bash
pip install -r requirements.txt      # or requirements-lock.txt for exact versions
streamlit run app.py
```
For structure prep, also install Open Babel (`apt install openbabel`, or a conda env
on Windows — see the guides below).

## Deploy it
Self-host on your own server with **[SELF_HOSTING.md](SELF_HOSTING.md)** (Ubuntu +
nginx + HTTPS + systemd), or use Streamlit Community Cloud (**[DEPLOY.md](DEPLOY.md)**).

## Verify it
```bash
python -m tests.test_core        # 23 tests — engine + stats + figures
# full suite: 66 tests across 9 files, all validated against ground truth
```

---

## Project layout
| File | Role |
|---|---|
| `neurofig_core.py` | Deterministic engine: load, infer, choose+run tests, draw, export |
| `doseresponse.py` / `standardcurve.py` / `cluster_stats.py` | Validated analysis methods |
| `gasteiger.py` / `pdbqt_tools.py` | Cheminformatics (charges, structure prep) |
| `labcalc.py` | Bench calculators |
| `journal_figures.py` | The 20-figure gallery (+ real-data versions) |
| `licensing.py` | Offline signed license keys (paywall) |
| `app.py` | Streamlit UI — calls the core only |
| `tests/` | 66 tests validating every component |
| `VALIDATION.md` | Full validation & reproducibility dossier |
| `SELF_HOSTING.md` | Personal-server deployment runbook |
| `NeuroFig_Blueprint.md` | Product/strategy blueprint |

---

## Scientific integrity
NeuroFig reports the test it used, never prints `p = 0`, refuses to extrapolate a
fabricated value (returns NaN out of range), and states its limitations openly in
[VALIDATION.md](VALIDATION.md). Uploaded data is processed in memory and not stored
after the session.
