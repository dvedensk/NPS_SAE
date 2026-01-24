## Model-level summary (ordered by MSE)

``` r
load("data/ACS_NPS_simulation_results.RData")
stopifnot(exists("summary_df"))

summary_df <- summary_df %>%
  arrange(MSE)

kable(
  summary_df,
  digits = c(5, 5, 5, 3, 3),
  caption = "Model metrics aggregated across simulations"
)
```

| model | MSE | MAB | Coverage | Int. Score |
|:---|---:|---:|---:|---:|
| bulm_uncond | 0.00083 | 0.02301 | 0.868 | 0.149 |
| bulm_beta | 0.00097 | 0.02529 | 0.836 | 0.159 |
| mrp-p | 0.00103 | 0.02628 | 0.826 | 0.163 |
| mrp-p-int | 0.00104 | 0.02645 | 0.833 | 0.166 |
| bulm_ps_only | 0.00214 | 0.03739 | 0.957 | 0.207 |
| IPW HT (uncond) | 0.00329 | 0.04888 | 0.673 | 0.375 |
| mrp-r-int-bootstrap | 0.00535 | 0.05686 | 0.932 | 0.367 |
| mrp-r-bootstrap | 0.00535 | 0.05703 | 0.929 | 0.367 |
| NPS Prior w/ Pseudolikelihood: Power Prior | 0.00583 | 0.05952 | 0.957 | 0.416 |
| NPS Prior w/ Raking: Power Prior | 0.00604 | 0.06073 | 0.950 | 0.409 |
| VSW | 0.00618 | 0.06270 | NA | NA |
| direst | 0.00873 | 0.07479 | 0.932 | 0.483 |
| IPW HT (beta) | 0.01933 | 0.13240 | 0.060 | 3.006 |

Model metrics aggregated across simulations

## Highlights

- **Top MSE:** bulm_uncond, bulm_beta, mrp-p, mrp-p-int (all ≲0.0011).
- **Coverage leaders:** bulm_ps_only (0.957), NPS prior PP
  (0.955–0.950), mrp-r variants (~0.93).
- **Interval Score (lower is better):** bulm_uncond/bulm_beta ≈0.15–0.16
  are best; IPW HT (beta) is worst (3.006).
- **Lowest performers (MSE):** IPW HT (beta) ≫ direst \> VSW.

## Notes

- VSW reports NA for interval metrics by design.
- NPS power-prior variants now produce intervals (pseudolikelihood and
  raking weight-covariate).
