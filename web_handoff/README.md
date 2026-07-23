# NeuroFig — Web Integration Handoff Package

**For the web developer building the scientific-services website.** This folder
contains everything you need to plan and incorporate NeuroFig into the site. Read
in this order:

| File | What it gives you |
|---|---|
| `01_product_brief.md` | What NeuroFig is, what it does, who it's for |
| `02_integration_guide.md` | **How to technically wire it into the website** (the important one) |
| `03_site_structure.md` | Recommended pages / information architecture |
| `04_brand_and_design.md` | Colours, type, tone — so the site matches the figures |
| `05_website_copy.md` | Ready-to-use copy (hero, tool cards, trust section, CTAs) |
| `06_business_model.md` | Freemium + licensing + the games funnel |
| `assets/` | Real showcase figure images produced by the tool |

## The 30-second version
NeuroFig is a **validated, no-code tool that turns a data upload into a journal-
ready figure with the correct statistics** — plus dose-response/IC50, standard
curves, structure prep (PDB→PDBQT), and 20 journal figure types. It is a self-
contained **Streamlit web app** (Python). The website's job is the marketing,
brand, SEO, pricing, and a **"Launch tool"** button that opens the app (hosted on
a subdomain). You do **not** need to rebuild the tool — just link to it and sell it.

## The one thing that makes this different
Everything is **deterministic and validated against ground truth** — no AI in the
compute path, nothing hallucinated. There is a public validation dossier
(`VALIDATION.md` in the main repo). **Lead with that trust** — it's the moat.
