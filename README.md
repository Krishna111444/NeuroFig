# NeuroFig — starter codebase

Turn raw experimental data into a publication-ready figure with the correct
statistics, in the browser, with no coding for the end user. This folder is a
**working MVP seed** you extend with Claude Code.

## What's here

| File | Role |
|------|------|
| `neurofig_core.py` | The deterministic engine: load data, infer columns, choose + run the right statistical test, render the figure. **No UI, no AI — testable on its own.** |
| `app.py` | The Streamlit web UI that wraps the engine: upload → confirm columns → figure + stats → download PNG/PDF/SVG. |
| `requirements.txt` | Python dependencies. |
| `_test_core.py` | (optional) a smoke test showing the engine works end-to-end. |

## Run it locally (5 minutes)

```bash
# 1. get the code onto your machine, then inside this folder:
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Your browser opens at `http://localhost:8501`. Tick "use a sample dataset" and
click **Generate figure** to see the whole loop with no data of your own.

## Deploy it free (so others can use it)

The fastest path is **Streamlit Community Cloud**: push this folder to a GitHub
repo, connect the repo at streamlit.io, and it hosts the app on a public URL for
free. Alternatives: Hugging Face Spaces, Render, Railway. None require a card to
start.

## How to work with Claude Code on this

Open Claude Code in this folder and extend it one function at a time. The code is
deliberately split so you can point Claude Code at a single piece. Good first asks:

- *"In `neurofig_core.py`, add a `run_nonparametric` option (Mann-Whitney for 2
  groups, Kruskal-Wallis + Dunn for 3+) and let `app.py` offer it when the Shapiro
  warning fires."*
- *"Add a second figure type to `neurofig_core.py`: a paired before/after plot with
  connecting lines, and a Wilcoxon signed-rank test."*
- *"In `app.py`, add journal style presets (Nature, Cell, eLife) that change font
  sizes, figure width in mm, and colour palette."*
- *"Add a watermark to the free PNG export and gate the clean PDF/SVG behind a
  simple payment check."*

Always ask Claude Code to **run `_test_core.py` after each change** so you catch
breakage immediately even though you don't read the code yourself.

## The one cost to remember

Running this app itself is free/cheap. The **AI layer** you add later (plain-English
→ figure) calls a model API that is billed per use — that is your only real running
cost, and your price per figure must stay above it. See the blueprint document for
the full monetization plan.
