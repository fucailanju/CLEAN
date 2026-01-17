# Updated by Jack on 01162026
# The test version when SBM is obvious yet sparse

library(CLEAN)
set.seed(1)

ids10 <- paste0("i", 1:10)
A10 <- matrix(0, 10, 10, dimnames = list(ids10, ids10))

# two dense blocks: 1-5 and 6-10
A10[1:5, 1:5] <- 1
diag(A10[1:5, 1:5]) <- 0
A10[6:10, 6:10] <- 1
diag(A10[6:10, 6:10]) <- 0

# very sparse between-block ties (optional, keep it tiny)
A10[1,6] <- A10[6,1] <- 1

b10 <- clean_blocks(A10, seed = 1)
b10$degenerate
table(b10$memberships)
b10$memberships