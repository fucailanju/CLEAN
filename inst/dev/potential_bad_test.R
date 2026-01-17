# Updated by Jack on 01162026
# Potential bad tests

library(CLEAN)

# A) Degenerate network
A0 <- matrix(0, 6, 6, dimnames = list(ids, ids))
res0 <- clean_blocks(A0)
res0$degenerate
res0$memberships

# B) Wrong shape should error
A_bad <- matrix(0, 3, 4)
try(clean_validate_network(A_bad))

# C) Missing time should error
A_missing <- list("2000" = A)
try(clean_blocks_by_time(df, "id", "year", A_missing))
