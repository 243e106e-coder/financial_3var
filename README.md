# 3-variable financial TVP-GVAR branch

This folder is a clean branch of the existing 14-economy project. It does **not** overwrite the legacy 5-variable / trade-weight pipeline under `8.12/`.

## Research specification in this branch

- Economies: `AU BR CA CH CN EA UK JP KR NO SG TR US ZA`
- Domestic financial variables: `r`, `de`, `deq`
  - `r`: short-term interest-rate level
  - `de`: REER log change
  - `deq`: equity return
- Financial network: `GCAP_financial_W_2017.csv`
- Weight concept: GCAP 2017 bilateral portfolio positions, residence basis, row-normalized within the 14-economy sample
- Global risk variable: global GPR (`LN_GPR_QMEAN`)
- Brent is retained as a global control when the existing Brent file is available
- Default common sample: `2002Q2-2022Q4`

## Why this is a separate branch

The existing code hard-codes five domestic variables and the legacy trade-weight matrix. The purpose here is to change only the two dimensions central to the new financial-market paper:

1. 5 variables -> 3 financial-market variables;
2. trade W -> portfolio-financial W.

The first run deliberately keeps `p=1`, `q=1` and the existing global controls. This makes the stability comparison interpretable.

## Files

- `20_build_financial_3var_input.R`: extracts the 3 domestic variables, validates the 2017 GCAP W, constructs `r*`, `de*`, `deq*`, and merges global controls.
- `21_financial_weight_preflight.R`: validates the 14x14 matrix and compares it with the old trade W when available.
- `22_financial_3var_stability.R`: estimates an OLS diagnostic GVAR and calculates local and global spectral radii, `G0` conditioning, and a financial-vs-trade network comparison.
- `GCAP_financial_W_2017.csv`: fixed 2017 financial weight matrix.

## Run order

```bash
Rscript financial_3var/20_build_financial_3var_input.R
Rscript financial_3var/21_financial_weight_preflight.R
Rscript financial_3var/22_financial_3var_stability.R
```

Or manually run the GitHub Actions workflow **Financial 3var preflight**.

## Decision rule after the first run

Review `financial_3var/results/stability/01_global_stability_comparison.csv` first.

- If the financial specification has `GlobalSpectralRadius < 1`, proceed to the Bayesian TVP-GVAR stage.
- If it remains >= 1, diagnose the country/variable persistence before changing lag order, identification, weights, and transformations simultaneously.

This first-stage code performs no ad-hoc interpolation of the macro data or financial W.
