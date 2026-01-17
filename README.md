# CLEAN

CLEAN provides a lightweight pipeline to proxy **latent homophily structure** (e.g., via stochastic block model recovery) in networked panel data.

The package focuses on **structural preprocessing**, not estimation:

-   validating adjacency matrices,

-   recovering latent blocks via SBM,

-   attaching time-specific block proxies to panel data.

CLEAN is designed to be transparent, modular, and safe for research workflows.

------------------------------------------------------------------------

## Installation (development version)

``` r
# install.packages("remotes")
remotes::install_github("fucailanju/CLEAN")
```

------------------------------------------------------------------------

## Example 1: Block recovery

``` r
library(CLEAN)

ids <- paste0("i", 1:6)

A <- matrix(0, 6, 6, dimnames = list(ids, ids))
A[1,2] <- A[2,1] <- 1
A[1,3] <- A[3,1] <- 1
A[4,5] <- A[5,4] <- 1
A[4,6] <- A[6,4] <- 1

res <- clean_blocks(A)
res$memberships
```

------------------------------------------------------------------------

## Example 2: Panel data with temporally-varying blocks

``` r
set.seed(1)

df <- data.frame(
  id   = rep(ids, times = 2),
  year = rep(c(2000, 2001), each = 6),
  y    = rnorm(12),
  Wy   = rnorm(12),
  x1   = rnorm(12)
)

A_by_time <- list(
  "2000" = A,
  "2001" = A
)

df2 <- clean_blocks_by_time(
  df = df,
  id = "id",
  time = "year",
  A_by_time = A_by_time,
  model = "bernoulli"
)

head(df2[, c("id", "year", "clean_block")])
```

------------------------------------------------------------------------

## Notes

-   Latent blocks are recovered separately for each time period.
-   Block labels are time-suffixed to avoid cross-period contamination.
-   If a network is degenerate (no edges or a single block), observations are assigned to "no_block" for that period.
