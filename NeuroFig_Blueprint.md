# NeuroFig — Innovator's Blueprint & Build Portfolio

*A pre-submission figure-and-statistics engine for researchers.*
Prepared for Krishna · July 2026

This is the plan an innovator would hand a builder. It contains the product
concept, the architecture, exactly how the code fits together, the build process
step by step (with the prompts to give Claude Code), the monetization model, and
the risks. It is paired with a working starter codebase (`app.py`,
`neurofig_core.py`) so you can go from reading to running today.

---

## 1. The idea in one paragraph

Researchers waste hours turning raw data into figures that meet journal standards,
and the tools force a bad trade: matplotlib and R are powerful but need coding;
GraphPad Prism is no-code but costs money and is rigid. **NeuroFig** collapses that
trade. A scientist uploads their data (or describes what they want in plain
English), and the tool runs the *correct* statistical test automatically and returns
a colourblind-safe, vector, journal-ready figure — no code, no guesswork about which
test to use. The statistics and drawing are done by **deterministic code** (so they
can't be hallucinated); an AI layer sits on top only to translate plain English into
the right analysis. Sold per-figure and by subscription.

**Why this and not the "review papers" idea:** figure generation is deterministic,
so it is trustworthy, demoable in seconds, and safe to sell. It is a sharp wedge you
can win alone. The review / stats-audit / grant modules become later tabs of the
same product once you have paying users.

---

## 2. Who pays, and why

The buyer is a PhD student, postdoc, or PI preparing a manuscript, poster, or
thesis. Their pain is acute and deadline-driven: a figure that a reviewer calls
"unclear" or "not to standard" costs a revision round. They already pay for Prism
(~institutional licences) or fight matplotlib. A tool that produces a correct,
beautiful, correctly-tested figure in 30 seconds is worth a few dollars every time.

Start in **one discipline you can reach** — neuroscience is a strong choice: it is
figure-heavy, statistics-heavy, has strong conventions (mean±SEM, individual points,
significance stars, ΔF/F traces), and has large online communities to launch into.

---

## 3. The architecture (how the pieces fit)

Think of the product as four layers stacked on top of each other. You build them
bottom-up; each layer is useful before the one above exists.

```
┌──────────────────────────────────────────────────────────────┐
│  LAYER 4 — AI INTERFACE (add last)                             │
│  Plain English -> intent. "box plot of my 3 groups, Nature     │
│  style, with significance stars" -> a structured request the   │
│  engine below understands. This is the only layer that calls   │
│  a paid model API.                                             │
├──────────────────────────────────────────────────────────────┤
│  LAYER 3 — WEB APP / ACCOUNTS / PAYMENTS  (app.py + additions) │
│  Upload, configure, preview, download. Login. Paywall on the   │
│  clean vector export. Usage metering.                          │
├──────────────────────────────────────────────────────────────┤
│  LAYER 2 — DETERMINISTIC ENGINE  (neurofig_core.py) ⭐ DONE    │
│  Load data -> infer columns -> pick + run the correct test ->  │
│  render a journal figure -> export PNG/PDF/SVG. Pure Python,   │
│  no AI, fully testable. THIS IS THE PRODUCT'S HEART.           │
├──────────────────────────────────────────────────────────────┤
│  LAYER 1 — RENDERING + STATS PRIMITIVES (borrowed, never built)│
│  matplotlib (drawing) + scipy/statsmodels (statistics). Free,  │
│  battle-tested. You stand on these; you do not reinvent them.  │
└──────────────────────────────────────────────────────────────┘
```

The starter codebase already implements Layers 1–2 and a basic Layer 3. Your job
with Claude Code is to deepen Layer 3 and add Layer 4.

### How the code is organised (and why)
- **`neurofig_core.py`** is the engine and has *no* Streamlit or AI code. That
  separation is deliberate: it can be tested in isolation (`_test_core.py`), reused
  behind any interface later (web, API, plugin), and extended one function at a time.
- **`app.py`** is only the interface. It calls the engine. Keeping UI and logic apart
  is the single most important habit for a non-coder working with Claude Code —
  you can point it at one file/function without risking the rest.

---

## 4. The build roadmap (milestones, not months)

Each milestone is shippable and testable. Do not start the next until the current
one runs.

**M0 — Run the seed (today).** Get `app.py` running locally, generate a figure from
the sample data. You now have a working core. *Deliverable: it runs on your machine.*

**M1 — Make the engine robust (week 1–2).** Handle messy real uploads: odd column
names, missing values, wide-format tables, non-normal data (offer non-parametric
tests when the normality warning fires). *Deliverable: it survives 10 real files
from labmates without crashing.*

**M2 — More figure types (week 2–3).** Add the forms neuroscience actually needs:
grouped/mean±SEM bar (done), box/violin, paired before-after with lines, and the
time-course/ΔF/F trace with SEM shading. Each is one function in the engine + one
option in the UI. *Deliverable: covers ~80% of figures in a typical paper.*

**M3 — Journal presets + polish (week 3–4).** Style presets (Nature/Cell/eLife
sizing and fonts), editable axis labels, figure-width-in-mm export. *Deliverable:
output that passes a journal's figure checklist.*

**M4 — Accounts, paywall, deploy (week 4–6).** Add login, meter usage, gate the
clean vector export behind payment, deploy to a public URL. *Deliverable: a stranger
can pay and download.*

**M5 — The AI layer (week 6–8).** Add the plain-English box: user types a request,
a model returns a structured spec, the engine executes it. This is the "wow" and the
only part with a per-use cost. *Deliverable: "type what you want, get the figure."*

**M6 — Second module (later).** Reuse the same upload + account plumbing to add the
statistical-consistency checker, then grants. The all-in-one vision, assembled
safely one paying module at a time.

---

## 5. The build *process* — how a non-coder drives Claude Code

You are the **pilot and product owner**, not the coder. The loop for every change:

1. **Describe the outcome, not the code.** e.g. *"When the data isn't normal, offer a
   Mann-Whitney test instead of the t-test, and note it in the caption."*
2. **Ask Claude Code to change one file/function** and to keep the UI/engine split.
3. **Ask it to run `_test_core.py`** (and add a new test for the new behaviour) so
   breakage surfaces immediately — you don't need to read code to trust a passing test.
4. **Click through the app yourself** and look at the figure. Your eyes are the final
   check the tests can't do.
5. **Commit** (Claude Code can do this) so you can always roll back.

Copy-paste starter prompts for Claude Code are in `README.md`. Give it this blueprint
and the codebase as context at the start of each session.

### What you still own (no code, but real)
Model API account (fuel for Layer 4), a payment account (Razorpay for India / Paddle
or Lemon Squeezy for global cards, both handle tax), a domain, a privacy policy
(researchers guard unpublished data — promise you don't train on it and delete after
processing), and customer replies. None require programming.

---

## 6. Monetization

- **Pay-per-figure:** ~$3–5 per clean vector export. Matches how the pain occurs
  (occasional, deadline-driven) and has the lowest friction for first-time buyers.
- **Subscription:** ~$12–20/month for unlimited figures — for labs and serial
  submitters. This is where recurring revenue lives.
- **Free tier:** watermarked, lower-res preview. Lets people try before paying and
  drives word-of-mouth. The paywall gates only the *clean* download.
- **Later:** institutional/lab licences (a whole department paying once is your best
  revenue), and per-module upsells (stats-audit, grants).

**Unit economics that make this work broke:** your only variable cost is the Layer-4
model API call (a few cents per figure) plus tiny hosting. Price one paid figure to
cover many times its API cost and you are profitable from customer #1. Layers 1–3
cost effectively nothing to run.

---

## 7. The moat (why a copycat can't just clone it)

The LLM is not your moat — anyone can call the same model. What you accumulate is:
the **library of journal presets and field conventions**, the **hard-won handling of
messy real data**, the **correct test-selection logic**, and a **reputation in one
discipline's community**. Every figure processed teaches you another edge case that
goes back into the engine. That curated, discipline-specific knowledge is what a
weekend cloner doesn't have.

---

## 8. Honest risks & how to counter them

- **Wrong statistics = lost trust.** Counter: deterministic tests, explicit
  assumption checks (normality/variance), and *refuse* to run a test that doesn't fit
  rather than producing something wrong. Always show which test was used.
- **"It's just a matplotlib wrapper."** True until you add the conventions, presets,
  test-selection, and messy-data handling. That knowledge layer *is* the product.
- **Data privacy.** Unpublished data is sensitive. State clearly: no training on user
  files, deletion after processing. Lean into this — it's a selling point.
- **Scope creep (the all-in-one trap).** Counter: ship the figure wedge first. Modules
  come only after paying users. This blueprint is sequenced to enforce that.
- **You depend on Claude Code to change things.** Counter: keep the tight test loop,
  keep UI/engine separate, and commit often so you can always roll back.

---

## 9. The one-line version

Build the deterministic figure-and-stats engine first (it's done — `neurofig_core.py`),
wrap it in a simple upload→figure→download web app, launch it to neuroscientists as
a $3-a-figure tool, then add the plain-English AI layer as the "wow." Your moat is
the accumulated journal-and-field knowledge, not the model — and your only real cost
is the AI calls, which your price covers from the first customer.
