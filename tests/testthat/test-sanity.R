# tests/testthat/test-sanity.R

test_that("clean_validate_network enforces square matrices", {
  A_bad <- matrix(0, 3, 4)
  expect_error(clean_validate_network(A_bad), "square")
})

test_that("clean_validate_network rejects NA values", {
  A <- matrix(0, 3, 3)
  A[1,2] <- NA
  expect_error(clean_validate_network(A), "NA")
})

test_that("clean_validate_network checks symmetry when directed = FALSE", {
  ids <- paste0("i", 1:4)
  A <- matrix(0, 4, 4, dimnames = list(ids, ids))
  A[1,2] <- 1
  # asymmetric on purpose
  expect_error(clean_validate_network(A, directed = FALSE), "symmetric")
})

test_that("clean_blocks flags degenerate all-zero networks", {
  ids <- paste0("i", 1:6)
  A0 <- matrix(0, 6, 6, dimnames = list(ids, ids))
  res0 <- clean_blocks(A0, seed = 1)
  expect_true(res0$degenerate)
  expect_true(all(is.na(res0$memberships)))
})

test_that("clean_blocks is deterministic given seed on a clear two-block network", {
  ids10 <- paste0("i", 1:10)
  A10 <- matrix(0, 10, 10, dimnames = list(ids10, ids10))

  # two blocks: 1-5 and 6-10
  A10[1:5, 1:5] <- 1
  diag(A10[1:5, 1:5]) <- 0
  A10[6:10, 6:10] <- 1
  diag(A10[6:10, 6:10]) <- 0

  # tiny between-block tie
  A10[1,6] <- A10[6,1] <- 1

  m1 <- clean_blocks(A10, seed = 123)$memberships
  m2 <- clean_blocks(A10, seed = 123)$memberships
  expect_identical(m1, m2)

  # sanity: should find > 1 block (labels may swap, but counts should be 5/5)
  expect_equal(sort(as.integer(table(m1))), c(5L, 5L))
})

test_that("clean_blocks_by_time errors if a time slice is missing in A_by_time", {
  ids <- paste0("i", 1:6)
  df <- data.frame(
    id   = rep(ids, times = 2),
    year = rep(c(2000, 2001), each = 6),
    y    = rnorm(12),
    Wy   = rnorm(12),
    x1   = rnorm(12)
  )

  A <- matrix(0, 6, 6, dimnames = list(ids, ids))
  A[1,2] <- A[2,1] <- 1

  A_missing <- list("2000" = A)  # no 2001
  expect_error(
    clean_blocks_by_time(df, "id", "year", A_missing),
    "No adjacency matrix found"
  )
})

test_that("clean_blocks_by_time attaches labels and respects id ordering per time", {
  ids <- paste0("i", 1:6)

  # panel has same ids each year
  df <- data.frame(
    id   = rep(ids, times = 2),
    year = rep(c(2000, 2001), each = 6),
    y    = rnorm(12),
    Wy   = rnorm(12),
    x1   = rnorm(12)
  )

  # year-specific networks
  A2000 <- matrix(0, 6, 6, dimnames = list(ids, ids))
  A2001 <- matrix(0, 6, 6, dimnames = list(ids, ids))

  # 2000: two blocks (1-3) and (4-6)
  A2000[1:3, 1:3] <- 1; diag(A2000[1:3, 1:3]) <- 0
  A2000[4:6, 4:6] <- 1; diag(A2000[4:6, 4:6]) <- 0

  # 2001: degenerate (no edges) -> should be "no_block"
  # (leave A2001 all zeros)

  A_by_time <- list("2000" = A2000, "2001" = A2001)

  df2 <- clean_blocks_by_time(
    df = df,
    id = "id",
    time = "year",
    A_by_time = A_by_time,
    out_col = "clean_block",
    seed = 1
  )

  # year 2001 should be all "no_block"
  expect_true(all(df2$clean_block[df2$year == 2001] == "no_block"))

  # year 2000 should NOT be all no_block (should find 2 communities)
  y2000 <- df2$clean_block[df2$year == 2000]
  expect_true(any(y2000 != "no_block"))

  # labels should be time-suffixed with "_2000" for the 2000 slice
  expect_true(all(grepl("_2000$", y2000)))
})
