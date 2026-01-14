# CLEAN

CLEAN provides a minimal pipeline to proxy latent homophily structure
(e.g., via stochastic block model recovery) and re-estimate contagion
effects in panel data models.

The package is intentionally lightweight and focuses on:
- validating adjacency matrices,
- recovering latent blocks via SBM,
- attaching time-specific block proxies to panel data,
- re-estimating contagion models with block fixed effects.

---

## Installation (development version)

```r
# install.packages("remotes")
remotes::install_github("fucailanju/CLEAN")
```

---

## Required inputs

To use CLEAN, users must provide:

- A **panel data frame** containing:
  - a unit identifier (e.g., country code),
  - a time variable (e.g., year),
  - an outcome variable,
  - a contagion regressor (e.g., `Wy`).

- `A_by_time`: a **named list of adjacency matrices**, where:
  - each element corresponds to a single time period,
  - list names match the time values in the panel data,
  - row and column names of each matrix match the unit identifiers.

CLEAN does **not** construct adjacency matrices for the user.

---

## Example

```r
library(CLEAN)

res <- clean(
  df = df,
  formula = y ~ Wy + x1 + x2,
  id = "ccode",
  time = "year",
  A_by_time = A_by_time,
  model = "bernoulli",
  fit_fun = lm
)

summary(res$model)
res$diagnostics
```

---

## Notes

- Latent blocks are recovered separately for each time period.
- Block labels are time-suffixed to avoid cross-period contamination.
- If a network is degenerate (no edges or a single block), observations
  are assigned to `"no_block"` for that period.
