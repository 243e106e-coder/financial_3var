# Financial 3-variable TVP-GVAR

Independent repository for the 14-economy financial-market branch of the TVP-GVAR project.

## Baseline first-stage specification

- Economies: `AU BR CA CH CN EA UK JP KR NO SG TR US ZA`
- Domestic variables:
  - `r`: short-term interest-rate level
  - `de`: REER log change
  - `deq`: equity return
- Financial network: GCAP 2017 bilateral portfolio positions, residence basis
- Global GPR: `LN_GPR_QMEAN`
- Brent retained as a global control in the first diagnostic
- Diagnostic lags: `p=1`, `q=1`
- Default common sample: `2002Q2–2022Q4`

## Repository structure

```text
R/
  20_build_financial_3var_input.R
  21_financial_weight_preflight.R
  22_financial_3var_stability.R
data/
  GCAP_financial_W_2017.csv
.github/workflows/
  financial_3var_preflight.yml
```

The workflow shallow-clones the legacy repository into `source_repo/` only to read the already-validated macro data, GPR, Brent and legacy trade weights. It does not modify the old repository.

## First decision

Run **Financial 3var preflight** in GitHub Actions and inspect:

`results/stability/01_global_stability_comparison.csv`

The key statistic is the global spectral radius under `GCAP_financial_2017`. Do not start the expensive Bayesian TVP stage until this diagnostic is reviewed.
