# 03 · Site Structure (information architecture)

Guiding principle: **the site must be legible in 5 seconds.** A researcher should
instantly see it's a serious, validated tools company. Keep the paid tools and the
games in clearly separate zones (see `06_business_model.md` for why).

## Top-level navigation
```
[ Logo ]   Tools    How we validate    Pricing    Learn / Play    Blog    [ Launch tool ]
```

## Pages

### 1. Home (landing)
- **Hero:** one-line value prop + a real figure image + primary CTA "Launch tool"
  and secondary "See how it works".
- **Tool grid:** cards for each tool (use `assets/` images + copy from `05_website_copy.md`).
- **Trust strip:** "66 validated tests · no fabrication · public validation dossier".
- **Social proof / use-cases** (add as you get users).
- Footer: privacy, contact, links.

### 2. Tools (the money pillar)
A page (or a page per tool) describing each service, each with a "Launch" button to
the app subdomain:
- Figures + statistics
- Dose-response (IC50/EC50)
- Standard curve
- Time-course (ΔF/F)
- 20 journal figure types (gallery)
- PDB → PDBQT
- Bench calculators

### 3. How we validate (the differentiator — do not skip)
Publish the trust story: deterministic engine, validated against ground truth, no
AI in the compute path, honest limitations. Summarize `VALIDATION.md`. **Few
competitors have this page — it converts skeptical scientists.**

### 4. Pricing
- Free (watermarked previews) · Pro (clean vector exports, ~monthly or yearly) ·
  optional Lab/Institution licence.
- Wired to a checkout provider (see `02_integration_guide.md`).

### 5. Learn / Play (the funnel — clearly separate)
Science games and educational content. Different visual framing so it never makes
the paid tools look like toys. This is top-of-funnel: traffic, virality, SEO, and
future researchers.

### 6. Blog / Resources (SEO)
"How to make a volcano plot", "Which statistical test should I use?", "qPCR standard
curve explained" — each ends with a CTA into the relevant tool. This is how
researchers will *find* you.

## Sitemap
```
/
/tools  (/tools/figures, /tools/dose-response, /tools/standard-curve, ...)
/validate
/pricing
/learn        (games + education, its own sub-section)
/blog
/privacy
app.yoursite.com   (the Streamlit tool — separate)
```
