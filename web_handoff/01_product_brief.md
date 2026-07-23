# 01 · Product Brief

## What it is
NeuroFig turns a researcher's raw data into a **publication-ready figure with the
correct statistics**, in the browser, with no coding. The statistics and drawing
are done by deterministic, tested code (no machine-learning, nothing hallucinated),
and every method is validated against ground truth.

## Who it's for
PhD students, postdocs, and PIs in the life sciences preparing manuscripts, posters,
and theses — people who currently fight matplotlib/R or pay for GraphPad Prism.
Primary wedge: **neuroscience & pharmacology**; the tools generalize to any wet lab.

## The pain it removes
Turning data into a figure that (a) uses the *right* statistical test and (b) meets
a journal's formatting standard is slow, error-prone, and a common cause of reviewer
revisions. NeuroFig does it correctly in ~30 seconds.

## Feature list (for tool cards on the site)
1. **Figures + statistics** — upload data, auto-selects the correct test (t-test /
   ANOVA / non-parametric), draws a colourblind-safe journal figure.
2. **Dose-response (IC50/EC50)** — 4-parameter logistic fit with confidence
   intervals (the Prism staple).
3. **Standard curve** — protein assays / qPCR (with efficiency) / ELISA, and it
   back-calculates unknown sample concentrations.
4. **Time-course (ΔF/F)** — calcium/photometry traces with a proper cluster-based
   permutation test.
5. **20 journal figure types** — volcano, Manhattan, raincloud, swimmer, Circos,
   Sankey, treemap, split violin, heatmaps, and more.
6. **PDB → PDBQT** — clean a protein/ligand structure and convert it for molecular
   docking (AutoDock/Vina).
7. **Bench calculators** — molarity, dilutions, DNA copy number, primer/protein
   properties, dosing, and more (all reference-checked).
8. **Journal presets** — Nature/Cell/eLife/Science fonts and **exact figure-width-
   in-mm** vector export (PDF/SVG).

## Proof points (use these as trust signals)
- 66 automated tests, all passing, validating every component against ground truth.
- Dose-response recovers known parameters to 4 significant figures.
- Own Gasteiger-charge algorithm matches Open Babel exactly on 12/13 molecules.
- Never prints an invalid "p = 0"; returns "no result" rather than a fabricated one.
- A public validation & reproducibility dossier (`VALIDATION.md`).

## Competitive framing
- **vs GraphPad Prism:** no licence cost, no install, no manual test-picking, and
  it's transparent about which test it ran and why.
- **vs matplotlib/R:** no coding required.
- **The moat:** validated correctness + accumulated field conventions, not the UI.
