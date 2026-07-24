# Sample datasets

Synthetic-but-realistic files to try every NeuroFig module. Upload the matching
file in the app, pick the columns shown, and you'll reproduce the validated result
listed below. All values are generated (no real subjects); safe to share.

| File | Module / tab | Columns to pick | Expected result |
|---|---|---|---|
| `1_group_comparison.csv` | Figures → Group comparison | group=`genotype`, value=`freezing_pct` | One-way ANOVA F(2,33)=49.1, p<0.001; WT≈Rescue, KO lower |
| `2_timecourse_dff.csv` | Figures → Time-course (ΔF/F) | time=`time_s`; A=`CSplus_*`; B=`CSminus_*` | Cluster test: significant difference ~0–2.7 s |
| `3_dose_response.csv` | Figures → Dose-response | conc=`concentration_M`, response=`percent_inhibition` | EC50 ≈ 1×10⁻⁶ M, Hill ≈ 1.0, R²≈0.996 |
| `4_standard_curve_qpcr.csv` | Figures → Standard curve (qPCR) | conc=`starting_copies`, signal=`Ct` | slope −3.32, **efficiency ≈ 100%**, R²≈0.999 |
| `5_volcano_deseq.csv` | Gallery → Volcano | FC=`log2FoldChange`, p=`pvalue`, label=`gene` | Volcano; 6 p=0 rows auto-floored (warning shown) |
| `6_raincloud.csv` | Gallery → Raincloud | group=`treatment`, value=`biomarker` | Raincloud of Sham/Mild/Severe |
| `7_correlation.csv` | Gallery → Correlation heatmap | select all 6 columns | Recovers IL6↔TNFα +0.73, IL10↔CRP −0.65, glucose↔insulin +0.46 |
| `8_protein_1CRN.pdb` | PDB → PDBQT (receptor) | — | Real protein (crambin); clean → rigid receptor PDBQT with Gasteiger charges |
| `9_ligand_aspirin.pdb` | PDB → PDBQT (ligand) | — | Real ligand; → flexible PDBQT, 3 rotatable bonds, charges assigned |
| `10_venn_sets.csv` | Gallery → Venn diagram | one column per set (member IDs) | 3-way overlap of RNAseq/Proteomics/ATACseq gene lists |
| `11_enrichment.csv` | Gallery → Enrichment dot-plot | `term`, `gene_ratio`, `count`, `p_adjust` | GO/KEGG-style enrichment dot-plot |
| `12_pca_matrix.csv` | Gallery → PCA / ordination | `sample`, `group`, `feature*` | Samples separate by group on PC1/PC2 |
| `13_composition.csv` | Gallery → Stacked composition | `sample` + one column per taxon | Relative-abundance stacked bars |

**Multiomics figures** (Venn, UpSet, enrichment, PCA, stacked composition) currently
render in the gallery from built-in sample data; these CSVs show the upload format.

Reference conversion outputs are included (`*.receptor.pdbqt`, `*.ligand.pdbqt`) so
you can compare what the tool should produce.

**Calculators** need no file — try e.g. molarity: 0.1 mol/L × 1 L × 58.44 g/mol → 5.844 g.

`8_protein_1CRN.pdb` is the experimental structure of crambin (RCSB PDB ID 1CRN),
included here only as a public-domain test input.
