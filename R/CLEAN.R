
#' Reset / align adjacency matrix dimnames safely
#'
#' @param A adjacency matrix
#' @param ids optional character vector of node ids in desired order
#' @param mode "auto" (default), "rename_only", or "reorder_if_possible"
#' @return matrix with consistent row/col names (and optionally reordered)
#' @export
clean_reset_dimnames <- function(A, ids = NULL,
                                 mode = c("auto", "rename_only", "reorder_if_possible")) {
  mode <- match.arg(mode)

  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A)) stop("`A` must be a matrix or coercible to matrix.")
  if (nrow(A) != ncol(A)) stop("`A` must be square.")

  # Drop weird attributes (haven) but keep dim/dimnames
  attributes(A) <- attributes(A)[intersect(names(attributes(A)), c("dim", "dimnames"))]

  n <- nrow(A)
  rn <- rownames(A)
  cn <- colnames(A)

  # Ensure dimnames exist and are identical if either exists
  if (is.null(rn) && is.null(cn)) {
    rn <- as.character(seq_len(n))
    cn <- rn
    rownames(A) <- rn
    colnames(A) <- cn
  } else if (is.null(rn) && !is.null(cn)) {
    rownames(A) <- cn
    rn <- cn
  } else if (!is.null(rn) && is.null(cn)) {
    colnames(A) <- rn
    cn <- rn
  } else if (!identical(rn, cn)) {
    stop("Row and column names exist but do not match (rn != cn).")
  }

  if (is.null(ids)) return(A)

  ids <- as.character(ids)
  if (length(ids) != n) stop("Length of `ids` must equal nrow(A).")

  # Decide behavior
  if (mode == "auto") {
    # If existing names look informative (not 1..n), and contain ids -> reorder.
    looks_index <- suppressWarnings(all(!is.na(as.integer(rn))) && identical(rn, as.character(seq_len(n))))
    can_reorder <- all(ids %in% rn)
    mode <- if (!looks_index && can_reorder) "reorder_if_possible" else "rename_only"
  }

  if (mode == "reorder_if_possible") {
    if (!all(ids %in% rn)) stop("Cannot reorder: some ids not found in rownames(A).")
    A <- A[ids, ids, drop = FALSE]
  } else {
    # rename only, keep current order
    rownames(A) <- ids
    colnames(A) <- ids
  }

  A
}



#' Validate adjacency matrix for CLEAN
#'
#' @param A adjacency matrix (matrix/data.frame coercible to matrix)
#' @param ids optional character vector of node ids (panel ordering)
#' @param directed logical; if FALSE, checks symmetry (within tolerance)
#' @param allow_isolates logical; if FALSE, errors if any node has degree 0
#' @param tol numeric tolerance for symmetry check
#' @param zero_diag logical; if TRUE, forces diag(A)=0
#' @param sym_action what to do if directed=FALSE and A not symmetric:
#'        "error" (default) or "symmetrize"
#' @return A validated numeric matrix
#' @export
clean_validate_network <- function(A,
                                   ids = NULL,
                                   directed = FALSE,
                                   allow_isolates = TRUE,
                                   tol = 1e-10,
                                   zero_diag = TRUE,
                                   sym_action = c("error", "symmetrize")) {

  sym_action <- match.arg(sym_action)

  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A)) stop("`A` must be a matrix or coercible to a matrix.")
  if (nrow(A) != ncol(A)) stop("`A` must be square (nrow == ncol).")

  suppressWarnings(storage.mode(A) <- "numeric")
  if (anyNA(A)) stop("`A` contains NA values; set missing ties explicitly to 0.")

  # Harmonize dimnames, and if ids provided, align consistently
  A <- clean_reset_dimnames(A, ids = ids, mode = "auto")

  # If ids provided, force correct order + names
  if (!is.null(ids)) {
    ids <- as.character(ids)

    # If A has names (it does after reset), enforce exact order
    if (length(ids) != nrow(A)) stop("Length of `ids` must match nrow(A).")

    # Re-label & reorder to ids order:
    # since A currently has "1..n", we treat current order as data order,
    # and simply rename to ids without permuting.
    # (If you want true lookup/reorder by existing names, add another mode later.)
    rownames(A) <- ids
    colnames(A) <- ids
  }

  if (zero_diag) diag(A) <- 0

  if (!directed) {
    max_asym <- max(abs(A - t(A)))
    if (max_asym > tol) {
      if (sym_action == "error") {
        stop("`A` is not symmetric but `directed = FALSE`.")
      } else {
        # symmetrize safely
        A <- (A + t(A)) / 2
        if (zero_diag) diag(A) <- 0
      }
    }
  }

  deg <- rowSums(A != 0) + colSums(A != 0)
  if (!allow_isolates && any(deg == 0)) stop("`A` contains isolate nodes.")

  A
}


#' Detect latent blocks using a simple LSM-style embedding + clustering
#'
#' This is a lightweight LSM implementation:
#' 1) convert adjacency to a symmetric dissimilarity matrix
#' 2) embed into d-dimensional space via cmdscale()
#' 3) cluster the embedded points via kmeans()
#'
#' @param A adjacency matrix (square)
#' @param d embedding dimension (default 2)
#' @param k number of clusters; if NULL, chosen automatically (elbow heuristic)
#' @param k_max maximum k to consider when choosing k automatically
#' @param directed logical; if TRUE, symmetrizes A for embedding
#' @param seed random seed
#' @return list with memberships and embedding
#' @export
clean_blocks_lsm <- function(A,
                             d = 2,
                             k = NULL,
                             k_max = 10,
                             directed = TRUE,
                             seed = 1) {

  A <- clean_validate_network(A, directed = directed)

  set.seed(seed)

  n <- nrow(A)
  if (n <= 1 || sum(A != 0, na.rm = TRUE) == 0) {
    mem <- rep(NA_integer_, n)
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, d = d, k = k, embedding = NULL, degenerate = TRUE))
  }

  # --- Step 1: create a symmetric similarity matrix for embedding ---
  S <- A
  if (directed) S <- (A + t(A)) / 2
  diag(S) <- 0

  # If S is constant / near-constant, embedding degenerates
  if (max(S, na.rm = TRUE) == min(S, na.rm = TRUE)) {
    mem <- rep(NA_integer_, n)
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, d = d, k = k, embedding = NULL, degenerate = TRUE))
  }

  # --- Step 2: convert to dissimilarity (higher tie => closer => smaller distance) ---
  # Scale to [0,1] then convert to distance
  S01 <- (S - min(S, na.rm = TRUE)) / (max(S, na.rm = TRUE) - min(S, na.rm = TRUE))
  D <- 1 - S01
  diag(D) <- 0

  # cmdscale expects a distance object or dist-like matrix
  dd <- as.dist(D)

  # Make sure d is feasible
  d_use <- min(as.integer(d), n - 1L)
  if (is.na(d_use) || d_use < 1L) d_use <- 1L

  emb <- tryCatch(
    stats::cmdscale(dd, k = d_use, eig = TRUE),
    error = function(e) NULL
  )

  if (is.null(emb) || is.null(emb$points)) {
    mem <- rep(NA_integer_, n)
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, d = d_use, k = k, embedding = NULL, degenerate = TRUE))
  }

  X <- emb$points
  if (is.vector(X)) X <- matrix(X, ncol = 1)

  # --- Step 3: choose k if not provided (simple elbow on tot.withinss) ---
  if (is.null(k)) {
    k_upper <- min(k_max, n - 1L)
    if (k_upper < 2L) {
      k <- 1L
    } else {
      wss <- rep(NA_real_, k_upper)
      for (kk in 1:k_upper) {
        fit <- stats::kmeans(X, centers = kk, nstart = 20)
        wss[kk] <- fit$tot.withinss
      }
      # elbow heuristic: maximize second difference (curvature)
      if (k_upper >= 3) {
        d1 <- diff(wss)
        d2 <- diff(d1)
        k <- which.max(-d2) + 1L
      } else {
        k <- 2L
      }
    }
  }

  k <- as.integer(k)
  if (is.na(k) || k < 1L) k <- 1L
  if (k > n) k <- n

  if (k == 1L) {
    mem <- rep(1L, n)
  } else {
    km <- stats::kmeans(X, centers = k, nstart = 50)
    mem <- km$cluster
  }

  if (!is.null(rownames(A))) names(mem) <- rownames(A)

  list(memberships = mem,
       d = d_use,
       k = k,
       embedding = X,
       eig = emb$eig,
       degenerate = FALSE)
}




#' Detect latent blocks using SBM or LSM
#'
#' @param A adjacency matrix (square)
#' @param method "sbm" (default) or "lsm"
#' @param model SBM model passed to sbm::estimateSimpleSBM (bernoulli/poisson/gaussian)
#' @param directed logical; TRUE allows directed adjacency (LSM embedding will symmetrize)
#' @param seed random seed
#' @param d LSM embedding dimension (default 2)
#' @param k LSM number of clusters; if NULL auto-selects
#' @param k_max LSM max k to consider if k is NULL
#' @param ... passed to sbm::estimateSimpleSBM when method="sbm"
#' @return list with memberships and fit details
#' @export
clean_blocks <- function(A,
                         method = c("sbm", "lsm"),
                         model = c("bernoulli","poisson","gaussian"),
                         directed = TRUE,
                         seed = 1,
                         d = 2,
                         k = NULL,
                         k_max = 10,
                         ...) {

  method <- match.arg(method)

  if (method == "lsm") {
    return(clean_blocks_lsm(A, d = d, k = k, k_max = k_max, directed = directed, seed = seed))
  }

  # --- SBM route ---
  model <- match.arg(model)
  A <- clean_validate_network(A, directed = directed)

  if (!requireNamespace("sbm", quietly = TRUE)) {
    stop("Package 'sbm' is required. Install it with install.packages('sbm').")
  }

  set.seed(seed)

  if (sum(A != 0, na.rm = TRUE) == 0 || nrow(A) <= 1) {
    mem <- rep(NA_integer_, nrow(A))
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, model = model, fit = NULL, degenerate = TRUE))
  }

  fit <- sbm::estimateSimpleSBM(netMat = A, model = model, ...)
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
#' @param method "sbm" (default) or "lsm"
#' @param model SBM model (default 'bernoulli') if method="sbm"
#' @param directed logical; TRUE for directed adjacency
#' @param d LSM embedding dimension (default 2) if method="lsm"
#' @param k LSM clusters; if NULL auto-selects
#' @param k_max LSM max k to consider if k is NULL
#' @param out_col name of output column to add
#' @param seed random seed
#' @param ... passed to sbm::estimateSimpleSBM when method="sbm"
#' @return df with out_col added
#' @export
clean_blocks_by_time <- function(df, id, time, A_by_time,
                                 method = c("sbm", "lsm"),
                                 model = "bernoulli",
                                 directed = TRUE,
                                 d = 2,
                                 k = NULL,
                                 k_max = 10,
                                 out_col = "clean_block",
                                 seed = 1,
                                 ...) {

  method <- match.arg(method)

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

    # Validate and align adjacency to this slice's id order
    A_t <- clean_validate_network(A_by_time[[tval]], ids = ids_t, directed = directed)

    # Degenerate slice => no_block
    if (sum(A_t != 0, na.rm = TRUE) == 0 || nrow(A_t) <= 1) {
      df[[out_col]][idx] <- "no_block"
      next
    }

    # One call only (SBM or LSM)
    res <- clean_blocks(
      A = A_t,
      method = method,
      model = model,
      directed = directed,
      seed = seed,
      d = d,
      k = k,
      k_max = k_max,
      ...
    )

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
#' @param formula regression formula for the outcome model (without block FE is fine)
#' @param id unit id column name (string)
#' @param time time column name (string)
#' @param A_by_time named list of adjacency matrices
#' @param model SBM model (default 'bernoulli')
#' @param block_col name of block column to create
#' @param add_block_fe if TRUE, adds factor(block_col) to the formula
#' @param fit_fun model fitting function; default stats::lm
#' @param ... passed to fit_fun
#' @return list with fitted model, augmented data, and diagnostics
#' @export
clean <- function(df, formula, id, time, A_by_time,
                  model = "bernoulli",
                  block_col = "clean_block",
                  add_block_fe = TRUE,
                  fit_fun = stats::lm,
                  ...) {

  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")
  if (!is.function(fit_fun)) stop("`fit_fun` must be a function.")

  # 1) add blocks
  df2 <- clean_blocks_by_time(
    df = df,
    id = id,
    time = time,
    A_by_time = A_by_time,
    model = model,
    out_col = block_col
  )

  # 2) optionally add block FE
  f2 <- formula
  if (add_block_fe) {
    # avoid double-adding if user already included it
    ftxt <- paste(deparse(formula), collapse = "")
    if (!grepl(block_col, ftxt, fixed = TRUE)) {
      f2 <- stats::as.formula(paste0(ftxt, " + factor(", block_col, ")"))
    }
  }

  # 3) fit model
  mod <- fit_fun(f2, data = df2, ...)

  # 4) minimal diagnostics
  diag <- list(
    n = nrow(df2),
    n_time = length(unique(as.character(df2[[time]]))),
    block_col = block_col,
    share_no_block = mean(df2[[block_col]] == "no_block", na.rm = TRUE)
  )

  list(model = mod, data = df2, formula = f2, diagnostics = diag)
}
