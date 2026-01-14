#' Validate adjacency matrix for CLEAN
#'
#' @param A adjacency matrix (matrix/data.frame coercible to matrix)
#' @param ids optional character vector of node ids (panel ordering)
#' @param directed logical; if FALSE, checks symmetry (within tolerance)
#' @param allow_isolates logical; if FALSE, errors if any node has degree 0
#' @param tol numeric tolerance for symmetry check
#' @return A validated numeric matrix
#' @export
clean_validate_network <- function(A,
                                   ids = NULL,
                                   directed = FALSE,
                                   allow_isolates = TRUE,
                                   tol = 1e-10) {

  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A)) stop("`A` must be a matrix or coercible to a matrix.")
  if (nrow(A) != ncol(A)) stop("`A` must be square (nrow == ncol).")

  suppressWarnings(storage.mode(A) <- "numeric")
  if (anyNA(A)) stop("`A` contains NA values; set missing ties explicitly to 0.")

  rn <- rownames(A)
  cn <- colnames(A)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    stop("Row and column names of `A` must match.")
  }

  if (!is.null(ids)) {
    ids <- as.character(ids)
    if (!is.null(rn)) {
      if (!all(ids %in% rn)) stop("Some `ids` not found in rownames(A).")
      A <- A[ids, ids, drop = FALSE]
    } else {
      if (length(ids) != nrow(A)) stop("Length of `ids` must match nrow(A).")
      rownames(A) <- ids
      colnames(A) <- ids
    }
  }

  if (!directed && max(abs(A - t(A))) > tol) {
    stop("`A` is not symmetric but `directed = FALSE`.")
  }

  deg <- rowSums(A != 0) + colSums(A != 0)
  if (!allow_isolates && any(deg == 0)) {
    stop("`A` contains isolate nodes.")
  }

  A
}


#' Detect latent blocks using SBM (minimal wrapper)
#'
#' @param A adjacency matrix (square)
#' @param model SBM model passed to sbm::estimateSimpleSBM (default 'bernoulli')
#' @param seed random seed
#' @return list with memberships and fit
#' @export
clean_blocks <- function(A, model = "bernoulli", seed = 1) {
  A <- clean_validate_network(A, directed = TRUE)

  if (!requireNamespace("sbm", quietly = TRUE)) {
    stop("Package 'sbm' is required. Install it with install.packages('sbm').")
  }

  set.seed(seed)

  # Degenerate network: no edges
  if (sum(A != 0, na.rm = TRUE) == 0 || nrow(A) <= 1) {
    mem <- rep(NA_integer_, nrow(A))
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, model = model, fit = NULL, degenerate = TRUE))
  }

  fit <- sbm::estimateSimpleSBM(netMat = A, model = model)
  mem <- fit$memberships
  if (!is.null(rownames(A))) names(mem) <- rownames(A)

  list(memberships = mem, model = model, fit = fit, degenerate = FALSE)
}


#' Detect blocks by time slice and attach to panel data
#'
#' @param df panel dataframe
#' @param id column name for unit id (string)
#' @param time column name for time (string)
#' @param A_by_time named list of adjacency matrices; names must match unique(df[[time]])
#' @param model SBM model (default 'bernoulli')
#' @param out_col name of output column to add
#' @param seed random seed
#' @return df with out_col added
#' @export
clean_blocks_by_time <- function(df, id, time, A_by_time,
                                 model = "bernoulli",
                                 out_col = "clean_block",
                                 seed = 1) {
  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!all(c(id, time) %in% names(df))) stop("`id` and/or `time` not found in df columns.")
  if (!is.list(A_by_time) || is.null(names(A_by_time))) {
    stop("`A_by_time` must be a *named* list of adjacency matrices.")
  }

  tt <- as.character(df[[time]])
  uu <- unique(tt)

  df[[out_col]] <- NA_character_

  for (tval in uu) {
    if (!tval %in% names(A_by_time)) {
      stop(sprintf("No adjacency matrix found for time '%s' in `A_by_time`.", tval))
    }

    idx <- which(tt == tval)
    ids_t <- as.character(df[[id]][idx])

    A_t <- clean_validate_network(A_by_time[[tval]], ids = ids_t, directed = TRUE)

    # If degenerate: assign no_block
    if (sum(A_t != 0, na.rm = TRUE) == 0 || nrow(A_t) <= 1) {
      df[[out_col]][idx] <- "no_block"
      next
    }

    res <- clean_blocks(A_t, model = model, seed = seed)
    mem <- res$memberships

    if (length(unique(mem[!is.na(mem)])) <= 1) {
      df[[out_col]][idx] <- "no_block"
    } else {
      df[[out_col]][idx] <- paste0("comm", mem, "_", tval)
    }
  }

  df
}

#' Apply CLEAN pipeline (minimal)
#'
#' @param df panel dataframe
#' @param id node id column
#' @param time time column
#' @param y outcome column
#' @param contagion contagion regressor column (already constructed by user)
#' @param A_by_time list of adjacency matrices
#' @param method block method
#' @return list with model + diagnostics
#' @export
clean <- function(df, id, time, y, contagion, A_by_time, method = "sbm") {
  stop("Not implemented yet")
}
