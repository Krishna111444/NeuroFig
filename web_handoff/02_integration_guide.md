# 02 · Integration Guide (technical)

**Read this first — it determines how you architect the site.**

NeuroFig is a **self-contained Streamlit app written in Python**. It is *not* a React
component you drop into the page, and it should not be rebuilt from scratch. The
validated Python engine is the asset; your job is to market it, sell it, and route
users to it.

## Recommended architecture: marketing site + app subdomain (standard SaaS)

```
yoursite.com            → the website YOU build (marketing, SEO, pricing, brand, games)
app.yoursite.com        → the NeuroFig Streamlit app (deployed separately)
```

- The **website** (your build — Next.js/Astro/WordPress/whatever) owns: landing
  pages, tool descriptions, the "How we validate" trust page, pricing/checkout,
  blog/SEO, and the games section.
- Each **"Launch tool"** button is just a link to `https://app.yoursite.com`.
- The Streamlit app is deployed on its own (a small VPS; the repo has a full runbook
  in `SELF_HOSTING.md`). It handles its own UI, uploads, figures, and paywall.

**Why this pattern:** Streamlit apps are not SEO-friendly and don't theme like a
marketing site — but they're perfect *as the tool*. Keeping the two separate is the
normal, clean way to do this (compare: linear.app the marketing site vs the actual
Linear app). You get a fast, beautiful, SEO-able site **and** a working tool without
fighting the framework.

## Integration options, ranked

1. **Link out to the app subdomain (recommended).** Simplest, cleanest UX, best
   performance. The button opens the app in a new tab or same tab.

2. **Embedded "mini demo" via iframe (optional, nice-to-have).** Streamlit supports
   an embed mode that hides its chrome:
   `https://app.yoursite.com/?embed=true`
   Put this in an `<iframe>` on a landing page as a live "try it right here" widget.
   Caveats: nested scrolling, limited mobile UX — good for a teaser, not the whole
   tool. Use link-out for the full experience.

3. **Native rebuild against an API (future, not now).** If you later want the tool
   fully inside your site's design system, wrap the Python engine (`neurofig_core.py`
   etc.) in a small FastAPI service and build a React UI against it. This is a
   substantial project — only pursue it once there's revenue justifying it. The
   engine is already cleanly separated from the UI, which makes this possible later.

## Payments & licensing (plan this on the website side)

The app already contains a working paywall — you don't build auth, but you build the
**pricing page + checkout**:

- Free tier: watermarked preview downloads (no login needed).
- Paid tier: a **license key** unlocks clean vector exports. Keys are signed and
  verified offline by the app — no user database required.
- **Flow to implement on the site:**
  1. Pricing page → checkout via **Lemon Squeezy / Paddle** (global cards, handles
     tax) or **Razorpay** (India/UPI). These give hosted checkout — you never touch
     card data.
  2. On successful payment, a license key is generated (`licensing.py` in the repo:
     `python licensing.py <email> <days>`) and emailed to the buyer. Start manual;
     automate later with the provider's "payment succeeded" webhook → a tiny endpoint
     that runs the key generator and emails it.
  3. Buyer pastes the key into the app's sidebar → unlocked.
- The app reads two secrets on the server: `NEUROFIG_LICENSE_SECRET` (signs/verifies
  keys) and `NEUROFIG_BUY_URL` (your checkout link, shown as the in-app Buy button).

## What the web developer needs to deliver
- The marketing website (pages in `03_site_structure.md`).
- A pricing page wired to a checkout provider.
- (Optional) an email step that sends the license key after purchase.
- DNS: an `app` subdomain A-record pointing to the app server.
- NOT: the tool itself, auth, or the figure/stat logic — those exist and are validated.

## Tech facts for planning
- App stack: Python 3.13, Streamlit 1.60, matplotlib/scipy/pandas/statsmodels;
  Open Babel (system package) for the docking-prep tab.
- Hosting: any small Linux VPS (2 GB RAM). Full deploy runbook: `SELF_HOSTING.md`.
- Data/privacy: uploads are processed in memory and not stored after the session —
  state this on the site; it's a selling point for unpublished data.
