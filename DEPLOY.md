# Deploying NeuroFig and taking payment

This is the runbook to get NeuroFig on a public URL where a stranger can pay and
download clean figures. Steps marked **(you)** require your own accounts/identity
and cannot be automated for you.

---

## 1. Run it locally first (sanity check)

```bash
pip install -r requirements.txt
streamlit run app.py
```
Open http://localhost:8501, tick the sample dataset, and generate a figure.

> If `pip install` fails with an SSL/certificate error, you are behind a network
> proxy intercepting HTTPS (seen on the SSSIHL campus network). Retry on a clean
> network or ask IT to whitelist `pypi.org`. The figure engine itself needs no network.

---

## 2. Put the code on GitHub **(you)**

```bash
git remote add origin https://github.com/<you>/neurofig.git
git push -u origin main
```
`.gitignore` already keeps your virtualenv, caches, and **secrets** out of the repo.

---

## 3. Deploy to Streamlit Community Cloud (free) **(you)**

1. Go to https://share.streamlit.io and sign in with GitHub.
2. "New app" → pick your repo, branch `main`, main file `app.py`.
3. Under **Advanced settings → Secrets**, paste (see `.streamlit/secrets.toml.example`):
   ```toml
   NEUROFIG_LICENSE_SECRET = "<a long random string>"
   NEUROFIG_BUY_URL = "<your checkout link>"
   ```
   Generate the secret once: `python -c "import secrets; print(secrets.token_urlsafe(32))"`.
4. Deploy. You now have a public `https://<app>.streamlit.app` URL.

Alternatives if you prefer: Hugging Face Spaces, Render, Railway — all work with the
same repo.

---

## 4. Take payment **(you)**

You need a checkout provider. Both below handle tax/invoicing and give you a hosted
payment page — **you never handle card numbers, and neither does this app.**

- **Global cards:** [Lemon Squeezy](https://www.lemonsqueezy.com) or Paddle
  (merchant-of-record: they remit tax for you).
- **India (INR/UPI):** [Razorpay](https://razorpay.com) Payment Links.

Create a product (e.g. "NeuroFig Pro — 1 year, $19"), copy its checkout URL, and set
it as `NEUROFIG_BUY_URL`.

---

## 5. Issue a license key after each sale **(you)**

When someone pays, mint them a key with your secret and email it:

```bash
NEUROFIG_LICENSE_SECRET="<same secret as in step 3>" \
  python licensing.py buyer@lab.edu 365
```
This prints a signed key valid for 365 days. The buyer pastes it into the app's
sidebar to unlock clean PDF/SVG exports. Keys are verified offline — no database.

> **Automate later:** point the provider's "payment succeeded" webhook at a small
> endpoint that runs `make_license(...)` and emails the key. Manual issuing is fine
> for your first customers.

---

## 6. Pricing (from the blueprint)

- Free: watermarked PNG preview (drives word-of-mouth).
- Pro: ~$12–20/mo or ~$19/yr for clean vector exports.
- Per-figure ($3–5) once you add metered checkout.

## What stays true about privacy
The app states, and you must honour: user data is processed only to make the figure,
is **not** used for training, and is **not** stored after the session. Researchers
guard unpublished data — this is a selling point, not fine print.
