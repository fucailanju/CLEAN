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

#' Detect latent blocks (minimal wrapper)
#'
#' @param A adjacency matrix
#' @param method block method, currently 'sbm'
#' @param K number of blocks or 'auto'
#' @param seed random seed
#' @return list with membership and metadata
#' @export
clean_blocks <- function(A, method = "sbm", K = "auto", seed = 1) {
  stop("Not implemented yet")
}

#' Detect blocks by time slice
#'
#' @param df panel dataframe
#' @param id node id column name
#' @param time time column name
#' @param A_by_time named list of adjacency matrices
#' @param method block method
#' @return df with a `clean_block` column
#' @export
clean_blocks_by_time <- function(df, id, time, A_by_time, method = "sbm") {
  stop("Not implemented yet")
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
