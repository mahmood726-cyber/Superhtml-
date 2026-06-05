# 786-MIII Pairwise (Synthesis)

A single-file, browser-based **pairwise meta-analysis** application. It computes
pooled effect sizes for binary (RR / OR / RD), continuous (MD / SMD), proportion
(logit), and survival (HR) data, with a random- or fixed-effects model and
several heterogeneity estimators. The interface includes a data editor, an
animated forest plot, funnel / Baujat / L'Abbe diagnostics, subgroup and
meta-regression views, sensitivity tools, a PRISMA flow builder, a clinical
(NNT/NNH) panel, a methods/results report generator, and a built-in validation
suite.

> The original file shipped with the title "786-MIII pair wise by Synthesis
> (Ultimate)". This README describes the tool by what it actually does; the word
> "Ultimate" is the historical file name, not a capability claim.

## Statistics

All pooling is performed **on the log scale** for ratio measures (RR, OR, HR)
and on the logit scale for proportions, then back-transformed for display.

- Effect sizes: log-RR, log-OR, RD, MD, Hedges' g (SMD, small-sample corrected),
  logit-proportion, log-HR (SE from the published CI width / 3.92).
- Pooling: inverse-variance and Mantel-Haenszel; fixed and random effects.
- Heterogeneity: Cochran's Q, I-squared, tau-squared via DerSimonian-Laird (DL),
  REML, Paule-Mandel (PM), and Sidik-Jonkman (SJ).
- Hartung-Knapp-Sidik-Jonkman (HKSJ) small-sample adjustment with the variance
  floor `max(1, q) * Var_RE` and a t-style critical value indexed by `k-1`.
- 95% confidence interval and (for random effects, k > 2) a prediction interval.
- Egger's regression test for funnel asymmetry.

The HKSJ critical value uses the approximation `t(df) ~= 1.96 + 2.37/df`
(`df = k-1`), not the exact Student-t quantile. The pooled-effect normal-tail
`p` uses the Zelen-Severo polynomial approximation. These are deliberate,
faithful design choices carried over from the original app, not bugs; they are
adequate for a teaching / exploratory tool and documented here for transparency.

## Offline / file layout

The app is fully offline — no CDN, no web fonts, no network calls.

- `index.html` — the application (UI + DOM/Plotly rendering inline). It is HTML
  despite the original file lacking an extension.
- `engine.js` — the pure statistical core (`Metric`, `Engine`), the single
  source of truth for all math. Loaded by `index.html` and required by the tests.
- `plotly.min.js` — vendored Plotly 2.24.1 (~3.5 MB), loaded locally.
- `tests.js` — pure Node test suite for `engine.js`.

## Running the tests

```
node tests.js
```

Exit code 0 means all checks passed. Every expected value in the suite is
hand-derived independently of the engine (see the comments in `tests.js`),
including a fully hand-worked two-study pooling example, k=1 passthrough,
two-identical-studies (tau-squared = I-squared = 0), empty guards, an HKSJ
floor + t-critical check, and the classic BCG-vaccine benchmark.

## Using the app

Open `index.html` in any modern browser (or visit the GitHub Pages URL). Choose a
data type, paste or enter data (or "Load Example"), then explore the Overview,
Forest, Diagnostics and other tabs.

## Fixes applied during revival (2026-06-05)

- **Offline vendoring:** replaced the Plotly CDN `<script>` with a locally
  vendored `plotly.min.js` (2.24.1); removed the Google Fonts `<link>` and
  preconnect tags. There are now zero external `http(s)` references.
- **Engine extraction:** moved the pure stats core (`Metric` + `Engine`) into
  `engine.js` verbatim and deleted the inline duplicate, so the math has a single
  source of truth shared by the app and the tests.
- **Tests added:** new `tests.js` (51 checks, all hand-derived) covering pooling,
  heterogeneity, HKSJ floor + t-critical, prediction-interval gating, and all
  four data types.
- **Pages-ready:** renamed the extension-less `superpairwise` file to
  `index.html`; added `.nojekyll`, `.gitignore`, and this README.

No statistical methodology was changed: the log-scale pooling, the HKSJ variance
floor, the t-indexed critical value, and the two-sided p-value were reviewed
against standard formulas and found correct/defensible, so the engine code was
copied verbatim rather than altered.

## License

MIT — see `LICENSE`.
