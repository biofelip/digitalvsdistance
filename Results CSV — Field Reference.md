# Results CSV — Field Reference

One row per simulated dataset. Every approach follows the same three-field
pattern: `P_c_X` (estimated average detection probability), `N_hat_X`
(estimated abundance), `bias_X` (estimation error). Reading `bias_X`:
**negative** = `P_c` underestimated → abundance overestimated; **positive**
= `P_c` overestimated → abundance underestimated.

## Identifiers

| Field | Meaning |
|---|---|
| `scenario` | Which of the 5 generation files this row came from (half-normal vs. beta `g_y`, low vs. high `r_y`, etc.) |
| `sim` | Simulation number within that scenario (1 to `nsim`) |

## Ground truth (known exactly, because this is simulated data)

| Field | Meaning |
|---|---|
| `n_detected` | Animals actually detected in this sim |
| `true_P_c` | **Exact** population-average detection probability — `mean(g)` over the *entire* simulated population (detected and undetected), not estimated from anything |
| `true_N_hat` | `n_detected / true_P_c` — sanity-check figure; should land close to the true population size (`n=10000`) |

## The seven fitted approaches

| # | Prefix | Model | What it represents |
|---|---|---|---|
| 1 | `plain` | Half-normal, distance only | Classical CDS — the textbook method. Assumes `g(0)=1` on the line, decay controlled purely by `sigma` |
| 2 | `cov` | Half-normal, `theta` on scale (log-link) | MCDS-style — same half-normal family, but `theta` stretches/shrinks the curve multiplicatively. Standard covariate approach |
| 3 | `strat` | Post-stratified plain half-normal | Cuts `theta` into `n_strata` quantile bins, fits approach 1 separately in each, pools via a Horvitz-Thompson-consistent weighted average — the classical "pooling robustness" workaround |
| 4 | `beta` | Beta density, distance only | Flexible shape (not constrained to monotonic decay from `y=0`), anchored at its own fitted mode rather than at the line |
| 5 | `beta_cov` | Beta density, `theta` on shape | Flexible-family analogue of approach 2 — `theta` modifies the beta's mean (`mu`) and concentration (`kappa`) via a log-link |
| 6 | `beta_theta` | Beta density fit to `theta` | Same flexible family, but fit directly to `theta` instead of distance — tests recoverability treating `theta` as the primary axis |
| 7 | `normal_theta` | Normal density fit to `theta` | Alternate flexible family for `theta`, for comparison against the beta version |

**`params_X`** (one text column per approach, format `"name=value; name=value"`) — the fitted parameter values, kept as a single compact string rather than one column per parameter so the file doesn't balloon in width:

- `params_plain`: `sigma_hat`
- `params_cov`: `beta_0`, `beta_1` (log-link intercept/slope on `sigma`)
- `params_beta`: `alpha`, `beta`
- `params_beta_cov`: `mu`, `gamma_0`, `gamma_1`
- `params_beta_theta`: `alpha`, `beta`
- `params_normal_theta`: `mu`, `sigma`
- *(approach 3 has none — it's a pooling procedure over approach-1 fits, not its own parametric model)*

**`n_strata_failed`** — how many of the `n_strata` quantile bins had too few animals or failed to fit, for approach 3 only.

## Baselines (not fitted models)

| Field | Meaning |
|---|---|
| `P_c_uniform` | Always exactly `1` — the strip-sampling assumption (perfect detection everywhere). `bias_uniform` directly equals `1 - true_P_c` |
| `N_hat_uniform` | Always equals `n_detected` — strip sampling never corrects for missed animals |
| `P_c_true_anchor0` | Uses the **true, exact** `g(y,theta)` surface (via a GAM smooth over the full population, no fitting error), rescaled by its value at `y=0` — the classical anchor, applied to the true shape |
| `P_c_true_anchormax` | Same true shape, but rescaled by its own **maximum** instead — the generalized anchor |
| `N_hat_true_anchor0` / `N_hat_true_anchormax` | Abundance from each of the above |

These two isolate bias coming purely from **which anchoring convention you pick**, with zero estimation error — since the shape is the true one, not a fitted one. Comparing `P_c_true_anchor0` against `P_c_true_anchormax` tells you how much the classical `g(0)=1` assumption alone costs you, independent of any model-fitting error.