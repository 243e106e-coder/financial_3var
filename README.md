# Financial 3-variable GVAR — clean 2017 financial-network baseline

This repository contains the clean first-stage pipeline for the 14-economy
financial branch of the TVP-GVAR project.

## Locked specification

- Economies: `AU BR CA CH CN EA UK JP KR NO SG TR US ZA`
- Domestic variables:
  - `r`: short-term interest-rate level
  - `de`: log change in the real effective exchange rate
  - `deq`: equity-market return
- Global variables: log GPR and Brent oil-price control
- Diagnostic lag order: `p=1`, `q=1`
- Main financial network: fixed 2017 GCAP nationality-restated portfolio
  exposure, `Restatement_Ex_Domestic`, equity plus bonds
- Robustness networks:
  - 2017 residency basis
  - 2017 nationality-restated equity only
  - 2017 nationality-restated bonds only

Rows of every weight matrix are investor economies and columns are issuer or
exposure-destination economies. The diagonal is zero and every row sums to one
within the 14-economy sample.

## Data integrity

No project script converts a missing bilateral position to zero, assumes
bilateral symmetry, or uses trade weights to fill financial positions.

The distributed GCAP data are not purely raw observations. The nationality
restatements are estimates, and the publisher documents that some residency
inputs can be interpolated or extrapolated. These publisher adjustments are
identified in `data/weights/cell_lineage_2017.csv` and
`data/weights/sources_and_limitations.csv`.

## Pipeline

```text
R/00_config.R
R/05_rebuild_weights_from_gcap.R     # independent source reconstruction
R/10_build_financial_3var_input.R    # common-sample domestic panel
R/20_weight_preflight.R              # four-matrix validation/comparison
R/30_stability_grid.R                # OLS GVAR stability under each W
R/40_estimation_gate.R               # READY/BLOCKED decision
```

Run **Financial 3var clean preflight** first. The expensive Bayesian TVP stage
must not start unless `results/gate/ESTIMATION_GATE.txt` reports `READY` for the
main network.

Run **Rebuild GCAP weights audit** whenever the committed weight matrices or
the reconstruction code changes. It downloads the public GCAP source,
reconstructs all four matrices, and requires an exact numerical match within
the declared tolerance.

## Main outputs

- `results/weights/01_weight_validation.csv`
- `results/weights/02_top_counterparts.csv`
- `results/weights/03_weight_concentration.csv`
- `results/weights/04_network_distance.csv`
- `results/stability/01_global_stability_comparison.csv`
- `results/stability/02_local_stability.csv`
- `results/gate/ESTIMATION_GATE.txt`

## Sources

- GCAP data: https://www.antoniocoppola.org/data
- Dataset documentation: https://globalcapitalallocation.s3.us-east-2.amazonaws.com/Restated_Bilateral_External_Portfolios_ReadMe.pdf
- QJE methodology: https://doi.org/10.1093/qje/qjab014
- IMF PIP/CPIS: https://data.imf.org/en/Data-Explorer?datasetUrn=IMF.STA:PIP(5.0.0)
- US Treasury TIC: https://home.treasury.gov/data/treasury-international-capital-tic-system
- BIS banking statistics are reserved for a separate bank-network robustness
  exercise and are not added to the portfolio matrices.
