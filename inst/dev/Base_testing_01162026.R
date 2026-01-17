# Updated by Jack on 01162026
# This is a testing .r file for CLEAN's sanity check and overall check 
# before formal publishing on GitHub

library(CLEAN)

set.seed(1)

ids <- paste0("i", 1:6)

A <- matrix(0, 6, 6, dimnames = list(ids, ids))
A[1,2] <- A[2,1] <- 1
A[1,3] <- A[3,1] <- 1
A[4,5] <- A[5,4] <- 1
A[4,6] <- A[6,4] <- 1

# basic: block detection
b <- clean_blocks(A)
b$degenerate
b$memberships

# panel example for blocks-by-time
df <- data.frame(
  id   = rep(ids, times = 2),
  year = rep(c(2000, 2001), each = 6),
  y    = rnorm(12),
  Wy   = rnorm(12),
  x1   = rnorm(12)
)

A_by_time <- list("2000" = A, "2001" = A)

df2 <- clean_blocks_by_time(
  df = df,
  id = "id",
  time = "year",
  A_by_time = A_by_time,
  model = "bernoulli",
  out_col = "clean_block",
  seed = 1
)

table(df2$clean_block, df2$year)
head(df2[order(df2$year, df2$id), c("id","year","clean_block")], 12)