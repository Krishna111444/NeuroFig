# 06 · Business Model & Strategy (context for planning)

This gives the web developer the *why* behind the structure, so decisions align with
the business.

## Revenue model — freemium
- **Free:** watermarked preview figures. Zero friction, drives word-of-mouth.
- **Pro (paid):** clean vector exports at exact journal widths, unlocked by a licence
  key. ~$12–20/month or ~$19–39/year (tune to market).
- **Lab / institution:** one licence for a whole group — the best revenue per sale.
- Variable cost is near zero (no per-figure AI cost, tiny hosting), so it's
  profitable from the first paying customer.

## How money flows (what the site implements)
1. Pricing page → hosted checkout (Lemon Squeezy/Paddle for global; Razorpay for
   India). The site never handles card data.
2. On payment → a signed licence key is generated and emailed (manual at first,
   webhook-automated later).
3. Buyer pastes the key into the tool → Pro unlocked. Keys verify offline; no user
   database needed.

## The two-audience strategy (why tools and games are separated)
| | Serious tools | Science games |
|---|---|---|
| Audience | Researchers (they pay) | Students, educators, curious public |
| Purpose | Correct answers; revenue | Reach, virality, education, SEO |
| Placement | Clean, credible "Tools" pillar | A distinct "Learn / Play" zone |

**The funnel logic:** games and educational blog content pull a large, sharable
audience → build brand and traffic → a fraction convert to tool users, and students
become tomorrow's researchers. This works **only if the two are visually separated** —
blended on one surface, the games make the paid tools look like toys and kill
credibility. Separate = funnel; blended = confusion.

## The moat (what the marketing should lean on)
Not the UI (anyone can copy a UI). The moat is:
1. **Validated correctness + transparency** — the "How we validate" page. Rare, and
   it converts skeptical scientists.
2. **Accumulated field conventions** — journal presets, correct test-selection, the
   messy-data handling — knowledge that compounds with every user.

## Positioning one-liner (for the whole business)
> Rigorous, validated scientific tools for researchers — the figures and analyses
> your reviewers won't bounce, with the receipts to prove they're correct.

## Founder focus note (honest strategy)
Early on, focus beats breadth. The validated tools are the hero and the revenue.
Games are a *supporting* funnel, not a co-equal pillar. Don't let the site become a
junk-drawer of everything — curate ruthlessly around the "trustworthy tools" promise.
