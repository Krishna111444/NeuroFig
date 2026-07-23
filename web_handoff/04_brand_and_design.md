# 04 · Brand & Design Tokens

Goal: the website must feel like the **same product** as the figures it produces.
The figures use the Okabe-Ito colourblind-safe palette — build the site palette from
the same colours so a figure dropped onto the page looks native.

## Colour palette (Okabe-Ito — the exact figure colours)
| Role | Hex | Notes |
|---|---|---|
| **Primary** | `#0072B2` | deep blue — buttons, links, brand |
| Accent / warning | `#D55E00` | vermillion — highlights, "up" states |
| Success | `#009E73` | green |
| Secondary | `#CC79A7` | reddish-purple |
| Highlight | `#E69F00` | orange |
| Light blue | `#56B4E9` | |
| Yellow | `#F0E442` | sparingly |
| Neutral | `#7F7F7F` | grey text/lines |
| Text | `#1A1A1A` | near-black |
| Background | `#FFFFFF` / `#F7F7F7` | white / off-white panels |

Support **light and dark** mode. Keep large areas neutral; use the palette for
accents and data, exactly like the figures.

## Typography
- Clean humanist sans-serif: **Inter**, or system `-apple-system, Segoe UI, Helvetica,
  Arial`. This mirrors the figures' sans-serif (Arial/Helvetica) so type feels
  consistent.
- Generous line-height, comfortable measure (~65–75 chars). Scientific audiences
  value clarity over flourish.

## Logo / identity direction
- Wordmark "NeuroFig" (or the parent brand for the services site) in the primary
  blue. A simple mark suggesting *data → figure* (e.g. points resolving into a curve
  or bar) works well. Keep it minimal and credible — this is a scientific-instruments
  brand, not a consumer app.

## Voice & tone
- **Rigorous, honest, no hype.** The brand promise is *correctness and transparency.*
- Say what it does plainly; state limitations openly; never overclaim.
- Words that fit: *validated, correct, reproducible, journal-ready, no fabrication,
  the right test.* Words to avoid: *magic, AI-powered* (the compute path is
  deliberately AI-free — don't imply otherwise), *revolutionary.*

## UI feel
- Calm, precise, lots of whitespace. Data-forward: show real figures, not stock
  photos. Fast. Accessible (WCAG AA; the palette is already colourblind-safe).

## Ready assets
`assets/` contains real tool output to use on the site:
- `showcase_group_bar.png` — bar chart + stats
- `showcase_dose_response.png` — IC50/EC50 sigmoid
- `showcase_volcano.png` — volcano plot
- `showcase_raincloud.png` — raincloud
- `showcase_trace.png` — ΔF/F trace with significance bar

These are exported by the actual tool, so they're the truest brand assets you have.
