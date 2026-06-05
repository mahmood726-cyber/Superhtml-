# E156 Protocol — 786-MIII Pairwise (Synthesis)

- **Project:** 786-MIII Pairwise (Synthesis) — offline pairwise meta-analysis app
- **Revived:** 2026-06-05
- **Type:** Single-file HTML tool (browser, offline) + pure JS engine + Node tests
- **Dashboard:** https://mahmood726-cyber.github.io/Superhtml-/

## What changed

- Vendored Plotly 2.24.1 locally (`plotly.min.js`) and removed Google Fonts; zero external references remain.
- Extracted the pure statistical core (`Metric` + `Engine`) into `engine.js` verbatim and deleted the inline duplicate (single source of truth).
- Added `tests.js` (51 hand-derived Node checks; all pass) covering pooling, heterogeneity, HKSJ, prediction intervals, and all four data types.
- Renamed the extension-less `superpairwise` file to `index.html`; added `.nojekyll`, `.gitignore`, README.
- Reviewed the math: log-scale pooling, HKSJ variance floor `max(1,q)*Var_RE`, and t-indexed critical value were confirmed correct/defensible; no methodology was altered.

## Body (E156 draft — CURRENT BODY)

How reliably can a single offline browser tool reproduce standard pairwise meta-analysis results across binary, continuous, proportion, and survival outcomes? The tool ingests user-entered study tables for relative-risk, odds-ratio, mean-difference, standardised-mean-difference, logit-proportion, and hazard-ratio data with no network access. Effect sizes are pooled on the log or logit scale using inverse-variance or Mantel-Haenszel weighting under fixed- or random-effects models, with DerSimonian-Laird, REML, Paule-Mandel, or Sidik-Jonkman heterogeneity estimators. The primary estimand is the pooled random-effects effect size with its 95% confidence interval, tau-squared, I-squared, and a prediction interval, back-transformed to the natural scale. A 51-check Node suite, hand-derived independently of the engine, confirms the two-study worked example, two-identical-studies homogeneity, k=1 passthrough, and the Hartung-Knapp variance floor. Results reproduce the classic BCG-vaccine benchmark within normal-approximation tolerance. The HKSJ critical value is a t-style approximation suitable for teaching, not exact-quantile inference.

SUBMITTED: [ ]
