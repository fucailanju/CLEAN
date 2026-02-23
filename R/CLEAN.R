

#' Reset / align adjacency matrix dimnames safely
#'
#' @param A a single adjacency matrix
#' @param ids optional character vector of node ids in desired order
#' @param mode "auto" (default), "rename_only", or "reorder_if_possible"
#' @param zero_diag logical; if TRUE set diag(A)=0
#' @return matrix with consistent row/col names (and optionally reordered)
#' @export
clean_reset_dimnames <- function(A, ids = NULL,
                                 mode = c("auto", "rename_only", "reorder_if_possible"),
                                 zero_diag = FALSE) {
  mode <- match.arg(mode)

  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A)) stop("`A` must be an n by n matrix or coercible to matrix.")
  if (nrow(A) != ncol(A)) stop("`A` must be square (n by n).")

  # Coerce to base matrix; drops most foreign attributes (e.g., haven-labelled)
  A <- as.matrix(A)
  suppressWarnings(storage.mode(A) <- "numeric")

  n <- nrow(A)
  rn <- rownames(A)
  cn <- colnames(A)

  # Ensure dimnames exist and are identical if either exists
  if (is.null(rn) && is.null(cn)) {
    rn <- as.character(seq_len(n))
    rownames(A) <- rn
    colnames(A) <- rn
  } else if (is.null(rn) && !is.null(cn)) {
    rownames(A) <- cn
    rn <- cn
  } else if (!is.null(rn) && is.null(cn)) {
    colnames(A) <- rn
    cn <- rn
  } else if (!identical(rn, cn)) {
    stop("Row and column names exist but do not match (rn != cn).")
  }

  if (zero_diag) diag(A) <- 0

  # If no ids, we’re done (just a cleaned-up matrix)
  if (is.null(ids)) return(A)

  ids <- as.character(ids)
  if (length(ids) != n) stop("Length of `ids` must equal nrow(A).")
  if (anyDuplicated(ids)) stop("`ids` contains duplicates; cannot safely align adjacency matrix.")

  # Decide behavior
  if (mode == "auto") {
    is_seq_index <- identical(rn, as.character(seq_len(n)))
    can_reorder  <- all(ids %in% rn)
    mode <- if (!is_seq_index && can_reorder) "reorder_if_possible" else "rename_only"
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


#' Reset dimnames for a named list of adjacency matrices
#'
#' @param A_by_time named list of adjacency matrices
#' @param ids_by_time optional named list of ids vectors (same names as A_by_time)
#' @param mode passed to clean_reset_dimnames()
#' @param zero_diag passed to clean_reset_dimnames()
#' @return named list of cleaned matrices
#' @export
clean_reset_dimnames_list <- function(A_by_time, ids_by_time = NULL,
                                      mode = "auto",
                                      zero_diag = FALSE) {
  if (!is.list(A_by_time) || is.null(names(A_by_time))) {
    stop("`A_by_time` must be a named list of matrices.")
  }

  if (!is.null(ids_by_time)) {
    if (!is.list(ids_by_time) || is.null(names(ids_by_time))) {
      stop("`ids_by_time` must be a named list if provided.")
    }
    missing_ids <- setdiff(names(A_by_time), names(ids_by_time))
    if (length(missing_ids) > 0) {
      stop("ids_by_time missing names: ", paste(missing_ids, collapse = ", "))
    }
  }

  out <- vector("list", length(A_by_time))
  names(out) <- names(A_by_time)

  for (nm in names(A_by_time)) {
    A_nm <- A_by_time[[nm]]
    ids  <- if (!is.null(ids_by_time)) ids_by_time[[nm]] else NULL

    # nicer error message early
    if (!is.null(ids)) {
      if (length(ids) != nrow(as.matrix(A_nm))) {
        stop("For time '", nm, "': length(ids_by_time[['", nm, "']]) must equal nrow(A_by_time[['", nm, "']]).")
      }
    }

    out[[nm]] <- clean_reset_dimnames(A_nm, ids = ids, mode = mode, zero_diag = zero_diag)
  }

  out
}



#' Validate adjacency matrix for CLEAN (single slice)
#'
#' @param A adjacency matrix (matrix/data.frame coercible to matrix)
#' @param ids optional character vector of node ids in desired order (panel ordering)
#' @param directed logical; if FALSE, enforce symmetry (or repair via sym_action)
#' @param allow_isolates logical; if FALSE, errors if any node has total degree 0
#' @param tol numeric tolerance (used for symmetry + integer checks)
#' @param zero_diag logical; if TRUE, forces diag(A)=0
#' @param sym_action what to do if directed=FALSE and A not symmetric:
#'        "error" or "symmetrize"
#' @param sym_rule for symmetrize: "avg" (continuous), "max"/"min" (count), "or"/"and" (binary)
#' @param value_type how to validate edge values:
#'        "auto" (no strict check), "binary", "count", "nonneg", "any"
#' @param allow_negative logical; if FALSE, stops if any negative entry
#' @return validated numeric matrix
#' @export
clean_validate_network <- function(A,
                                   ids = NULL,
                                   directed = FALSE,
                                   allow_isolates = TRUE,
                                   tol = 1e-10,
                                   zero_diag = TRUE,
                                   sym_action = c("error", "symmetrize"),
                                   sym_rule = c("auto", "avg", "max", "min", "or", "and"),
                                   value_type = c("auto", "binary", "count", "nonneg", "any"),
                                   allow_negative = FALSE) {

  sym_action <- match.arg(sym_action)
  sym_rule   <- match.arg(sym_rule)
  value_type <- match.arg(value_type)

  if (is.data.frame(A)) A <- as.matrix(A)
  if (!is.matrix(A)) stop("`A` must be a matrix or coercible to a matrix.")
  if (nrow(A) != ncol(A)) stop("`A` must be square (nrow == ncol).")

  suppressWarnings(storage.mode(A) <- "numeric")
  if (anyNA(A)) stop("`A` contains NA values; set missing ties explicitly to 0.")
  if (!all(is.finite(A))) stop("`A` contains Inf/-Inf; adjacency must be finite.")

  # Harmonize dimnames + optional reorder/rename
  A <- clean_reset_dimnames(A, ids = ids, mode = "auto")

  if (zero_diag) diag(A) <- 0

  if (!allow_negative && any(A < 0, na.rm = TRUE)) {
    stop("`A` contains negative values. If intended, set allow_negative = TRUE.")
  }

  # ---- value-type validation (tolerant) ----
  if (value_type == "binary") {
    if (any(abs(A - round(A)) > tol, na.rm = TRUE)) {
      stop("`A` must be binary (0/1) for value_type='binary' (non-integer values detected).")
    }
    Ar <- round(A)
    if (any(!(Ar %in% c(0, 1)))) stop("`A` must be binary (0/1) for value_type='binary'.")
    A <- Ar
    if (zero_diag) diag(A) <- 0
  } else if (value_type == "count") {
    if (any(A < 0, na.rm = TRUE) || any(abs(A - round(A)) > tol, na.rm = TRUE)) {
      stop("`A` must be nonnegative integers for value_type='count'.")
    }
    A <- round(A)
    if (zero_diag) diag(A) <- 0
  } else if (value_type == "nonneg") {
    if (any(A < 0, na.rm = TRUE)) stop("`A` must be nonnegative for value_type='nonneg'.")
  }

  # ---- symmetry check / repair when undirected ----
  if (!directed) {
    D <- A - t(A)
    diag(D) <- 0
    max_asym <- max(abs(D))

    if (max_asym > tol) {
      if (sym_action == "error") stop("`A` is not symmetric but `directed = FALSE`.")

      rule <- sym_rule
      if (rule == "auto") {
        rule <- if (value_type == "binary") "or"
        else if (value_type == "count") "max"
        else "avg"
      }

      if (value_type == "count" && rule == "avg") {
        stop("sym_rule='avg' is incompatible with value_type='count' (would create non-integers). Use 'max' or 'min'.")
      }

      if (rule == "avg") {
        A <- (A + t(A)) / 2
      } else if (rule == "max") {
        A <- pmax(A, t(A))
      } else if (rule == "min") {
        A <- pmin(A, t(A))
      } else if (rule == "or") {
        A <- 1 * ((A != 0) | (t(A) != 0))
      } else if (rule == "and") {
        A <- 1 * ((A != 0) & (t(A) != 0))
      } else stop("Unknown sym_rule.")

      if (zero_diag) diag(A) <- 0
    }
  }

  deg <- rowSums(A != 0) + colSums(A != 0)
  if (!allow_isolates && any(deg == 0)) stop("`A` contains isolate nodes.")

  A
}



#' Validate adjacency matrices for CLEAN (multiple time slices)
#'
#' Validates a *named list* of adjacency matrices (e.g., yearly/monthly slices).
#' For each slice name `nm`, this calls [clean_validate_network()] on `A_by_time[[nm]]`,
#' optionally aligning node labels/order using `ids_by_time[[nm]]`.
#'
#' @param A_by_time Named list of adjacency matrices. Each element must be a square matrix
#'   (or coercible to matrix). Names should match the time values used in the panel (e.g., "1951").
#' @param ids_by_time Optional named list of character vectors. If provided, it must contain
#'   the same names as `A_by_time`. Each `ids_by_time[[nm]]` gives the node IDs for that slice
#'   (panel ordering) and is passed as `ids=` into [clean_validate_network()].
#' @param ... Additional arguments forwarded to [clean_validate_network()], e.g.
#'   `directed`, `value_type`, `sym_action`, `sym_rule`, `zero_diag`, `allow_isolates`.
#'
#' @return A named list of validated adjacency matrices, with the same names as `A_by_time`.
#'
#' @examples
#' # Suppose A_by_time is a named list of yearly adjacency matrices:
#' # A_by_time <- list("2001" = A2001, "2002" = A2002)
#' # ids_by_time <- list("2001" = ids2001, "2002" = ids2002)
#' # A_ok <- clean_validate_network_list(A_by_time, ids_by_time,
#' #   directed = TRUE, value_type = "binary", zero_diag = TRUE)
#'
#' @export
clean_validate_network_list <- function(A_by_time, ids_by_time = NULL, ...) {
  if (!is.list(A_by_time) || is.null(names(A_by_time))) {
    stop("`A_by_time` must be a *named* list of adjacency matrices.")
  }

  if (!is.null(ids_by_time)) {
    if (!is.list(ids_by_time) || is.null(names(ids_by_time))) {
      stop("`ids_by_time` must be a *named* list when provided.")
    }
    missing <- setdiff(names(A_by_time), names(ids_by_time))
    if (length(missing) > 0) {
      stop("`ids_by_time` is missing keys for: ", paste(missing, collapse = ", "))
    }
  }

  out <- vector("list", length(A_by_time))
  names(out) <- names(A_by_time)

  for (nm in names(A_by_time)) {
    ids <- if (!is.null(ids_by_time)) ids_by_time[[nm]] else NULL

    out[[nm]] <- tryCatch(
      clean_validate_network(A_by_time[[nm]], ids = ids, ...),
      error = function(e) {
        stop(sprintf("Validation failed for slice '%s': %s", nm, e$message), call. = FALSE)
      }
    )
  }

  out
}


#' Latent Space Model (LSM) using latentnet::ergmm
#'
#' Fits an LSM and returns continuous latent positions (n x d).
#'
#' @param A adjacency matrix (square)
#' @param d latent dimension (default 2)
#' @param directed logical; whether network is directed
#' @param seed integer seed for reproducibility
#' @param ids optional character vector of node ids (panel ordering)
#' @param ergmm_control list passed to latentnet::control.ergmm()
#' @param ... passed to latentnet::ergmm()
#' @return list(positions = matrix n x d, fit = ergmm object, degenerate = logical)
#' @export
clean_lsm <- function(A,
                      d = 2,
                      directed = TRUE,
                      seed = 1,
                      ids = NULL,
                      ergmm_control = list(),
                      ...) {

  if (!requireNamespace("network", quietly = TRUE)) {
    stop("Package 'network' is required. Install it with install.packages('network').")
  }
  if (!requireNamespace("latentnet", quietly = TRUE)) {
    stop("Package 'latentnet' is required. Install it with install.packages('latentnet').")
  }

  # Validate + align names/order
  A <- clean_validate_network(A, ids = ids, directed = directed)

  n <- nrow(A)
  if (n <= 1 || sum(A != 0, na.rm = TRUE) == 0) {
    Z <- matrix(NA_real_, nrow = n, ncol = as.integer(d))
    if (!is.null(rownames(A))) rownames(Z) <- rownames(A)
    colnames(Z) <- paste0("lsm", seq_len(ncol(Z)))
    return(list(positions = Z, fit = NULL, degenerate = TRUE))
  }

  set.seed(seed)

  net <- network::network(
    A,
    directed = directed,
    matrix.type = "adjacency",
    loops = FALSE,
    ignore.eval = TRUE
  )

  # Sensible defaults; users can override via ergmm_control
  ctrl_default <- do.call(latentnet::control.ergmm, c(list(
    MCMCsamplesize = 2000,
    MCMCinterval   = 1024,
    burnin         = 2000
  ), ergmm_control))

  fit <- latentnet::ergmm(
    net ~ latentnet::euclidean(d = as.integer(d)),
    control = ctrl_default,
    ...
  )

  # ---- Extract positions Z robustly (latentnet versions differ) ----
  Z <- NULL

  # Common: sample$Z stores draws; take first draw (or last, if you prefer)
  Z <- tryCatch(fit$sample$Z[[1]], error = function(e) NULL)

  # Alternative slots seen in some versions
  if (is.null(Z)) Z <- tryCatch(fit$mkl$Z, error = function(e) NULL)

  if (is.null(Z)) {
    stop("Could not extract latent positions from ergmm fit (latentnet version mismatch).")
  }

  if (is.vector(Z)) Z <- matrix(Z, ncol = 1)
  if (nrow(Z) != n && ncol(Z) == n) Z <- t(Z)

  if (!is.null(rownames(A))) rownames(Z) <- rownames(A)
  colnames(Z) <- paste0("lsm", seq_len(ncol(Z)))

  list(positions = Z, fit = fit, degenerate = FALSE)
}



#' LSM by time slice: attach continuous latent coordinates to panel df
#'
#' @param df panel data.frame
#' @param id unit id column name (string)
#' @param time time column name (string)
#' @param A_by_time named list of adjacency matrices
#' @param d latent dimension (default 2)
#' @param directed logical
#' @param seed integer seed
#' @param out_prefix prefix for output columns (e.g., "clean_lsm_" -> clean_lsm_1, clean_lsm_2, ...)
#' @param ergmm_control list passed to control.ergmm()
#' @param on_fail what to do if a slice fails: "warn" (default) or "stop"
#' @param ... passed to clean_lsm() / ergmm()
#' @return df with new columns added
#' @export
clean_lsm_by_time <- function(df, id, time, A_by_time,
                              d = 2,
                              directed = TRUE,
                              seed = 1,
                              out_prefix = "lsm_",
                              ergmm_control = list(),
                              on_fail = c("warn", "stop"),
                              ...) {

  on_fail <- match.arg(on_fail)

  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!all(c(id, time) %in% names(df))) stop("`id` and/or `time` not found in df.")
  if (!is.list(A_by_time) || is.null(names(A_by_time))) {
    stop("`A_by_time` must be a *named* list of adjacency matrices.")
  }

  tt <- as.character(df[[time]])
  uu <- unique(tt)

  d <- as.integer(d)
  if (is.na(d) || d < 1L) stop("`d` must be a positive integer.")

  # Pre-create columns
  for (j in seq_len(d)) df[[paste0(out_prefix, j)]] <- NA_real_

  for (tval in uu) {

    if (!tval %in% names(A_by_time)) {
      msg <- sprintf("No adjacency matrix found for time '%s' in `A_by_time`.", tval)
      if (on_fail == "stop") stop(msg) else { warning(msg, call. = FALSE); next }
    }

    idx <- which(tt == tval)
    ids_t <- as.character(df[[id]][idx])

    # Try fit slice
    res <- tryCatch(
      clean_lsm(
        A = A_by_time[[tval]],
        d = d,
        directed = directed,
        seed = seed,
        ids = ids_t,
        ergmm_control = ergmm_control,
        ...
      ),
      error = function(e) e
    )

    if (inherits(res, "error")) {
      msg <- sprintf("LSM failed for slice '%s': %s", tval, res$message)
      if (on_fail == "stop") stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); next }
    }

    Z <- res$positions  # expected n x d

    # Degenerate / missing output => leave as NA and continue
    if (is.null(Z) || !is.matrix(Z) || nrow(Z) == 0L) {
      msg <- sprintf("LSM returned no positions for slice '%s' (leaving NA).", tval)
      if (on_fail == "stop") stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); next }
    }

    # Alignment safeguard: reindex by rownames(Z) if present
    if (!is.null(rownames(Z))) {
      missing_ids <- setdiff(ids_t, rownames(Z))
      if (length(missing_ids) > 0) {
        msg <- sprintf("Slice '%s': some ids missing from positions(): %s",
                       tval, paste(head(missing_ids, 10), collapse = ", "))
        if (on_fail == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
      }
      Z <- Z[ids_t, , drop = FALSE]  # will introduce NA rows if ids missing
    } else {
      # If no rownames, require exact length match to avoid silent misalignment
      if (nrow(Z) != length(ids_t)) {
        msg <- sprintf("Slice '%s': positions() has %d rows but ids_t has %d; cannot safely align.",
                       tval, nrow(Z), length(ids_t))
        if (on_fail == "stop") stop(msg, call. = FALSE) else { warning(msg, call. = FALSE); next }
      }
    }

    # Write back to df in idx order (ids_t order)
    for (j in seq_len(d)) {
      df[[paste0(out_prefix, j)]][idx] <- Z[, j]
    }
  }

  df
}



#' Detect latent blocks using SBM (single slice)
#'
#' @param A adjacency matrix (square)
#' @param model SBM likelihood model: "bernoulli", "poisson", or "gaussian"
#' @param directed logical; whether A should be treated as directed
#' @param seed integer; random seed for estimation
#' @param value_type how to validate edge values; default chosen from model if "auto"
#' @param sym_action what to do if directed=FALSE and A is not symmetric
#' @param sym_rule symmetrization rule when sym_action="symmetrize"
#' @param allow_isolates whether isolates are allowed
#' @param zero_diag force diag(A)=0
#' @param ... passed to sbm::estimateSimpleSBM
#' @return list(memberships, model, fit, degenerate)
#' @export
clean_sbm <- function(A,
                      model = c("bernoulli", "poisson", "gaussian"),
                      directed = TRUE,
                      seed = 1,
                      value_type = c("auto", "binary", "count", "nonneg", "any"),
                      sym_action = c("error", "symmetrize"),
                      sym_rule   = c("auto", "avg", "max", "min", "or", "and"),
                      allow_isolates = TRUE,
                      zero_diag = TRUE,
                      ...) {

  model      <- match.arg(model)
  value_type <- match.arg(value_type)
  sym_action <- match.arg(sym_action)
  sym_rule   <- match.arg(sym_rule)

  if (!requireNamespace("sbm", quietly = TRUE)) {
    stop("Package 'sbm' is required. Install it with install.packages('sbm').")
  }

  # Choose sensible defaults based on SBM model if user leaves auto
  if (value_type == "auto") {
    value_type <- switch(model,
                         bernoulli = "binary",
                         poisson   = "count",
                         gaussian  = "any"
    )
  }

  # Validate + (optionally) symmetrize if directed=FALSE
  A <- clean_validate_network(
    A = A,
    ids = NULL,
    directed = directed,
    allow_isolates = allow_isolates,
    zero_diag = zero_diag,
    value_type = value_type,
    sym_action = sym_action,
    sym_rule = sym_rule
  )

  # Degenerate: no edges or 1 node
  if (nrow(A) <= 1L || sum(A != 0, na.rm = TRUE) == 0L) {
    mem <- rep(NA_integer_, nrow(A))
    if (!is.null(rownames(A))) names(mem) <- rownames(A)
    return(list(memberships = mem, model = model, fit = NULL, degenerate = TRUE))
  }

  set.seed(seed)
  fit <- sbm::estimateSimpleSBM(netMat = A, model = model, ...)
  mem <- fit$memberships
  if (!is.null(rownames(A))) names(mem) <- rownames(A)

  list(memberships = mem, model = model, fit = fit, degenerate = FALSE)
}



#' Detect SBM blocks by time slice and attach to panel data
#'
#' @param df panel dataframe
#' @param id unit id column name
#' @param time time column name
#' @param A_by_time named list of adjacency matrices (names match unique(df[[time]]))
#' @param model SBM model: "bernoulli", "poisson", "gaussian"
#' @param directed logical; whether to treat each slice as directed
#' @param seed integer; base seed
#' @param seed_by_time logical; if TRUE, vary seed deterministically by time slice
#' @param out_col output column name
#' @param ... passed to sbm::estimateSimpleSBM
#' @return df with out_col added
#' @export
clean_sbm_by_time <- function(df, id, time, A_by_time,
                              model = c("bernoulli","poisson","gaussian"),
                              directed = TRUE,
                              seed = 1,
                              seed_by_time = FALSE,
                              out_col = "clean_block",
                              ...) {

  model <- match.arg(model)

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

    idx   <- which(tt == tval)
    ids_t <- as.character(df[[id]][idx])

    # Validate + align adjacency to df ordering for this time slice
    A_t <- clean_validate_network(
      A = A_by_time[[tval]],
      ids = ids_t,
      directed = directed,
      # value_type defaults will be chosen inside clean_sbm based on model,
      # but we want clean_validate_network() to be permissive here;
      # it will still enforce finiteness/NA/diag/names/order.
      value_type = "auto"
    )

    # Degenerate -> no_block
    if (nrow(A_t) <= 1L || sum(A_t != 0, na.rm = TRUE) == 0L) {
      df[[out_col]][idx] <- "no_block"
      next
    }

    this_seed <- if (seed_by_time) seed + as.integer(factor(tval)) else seed

    res <- clean_sbm(
      A = A_t,
      model = model,
      directed = directed,
      seed = this_seed,
      ...
    )

    mem <- res$memberships

    if (length(unique(mem[!is.na(mem)])) <= 1L) {
      df[[out_col]][idx] <- "no_block"
    } else {
      # IMPORTANT: memberships are slice-internal labels; that's OK.
      # They are only used for FE / controls within slice.
      df[[out_col]][idx] <- paste0("comm", mem, "_", tval)
    }
  }

  df
}



#' Apply CLEAN pipeline (SBM blocks or LSM coordinates)
#'
#' CLEAN = detect latent homophily from network(s), augment df, refit outcome model.
#'
#' @param df panel dataframe
#' @param formula outcome model formula (without homophily terms is fine)
#' @param id unit id column name (string)
#' @param time time column name (string)
#' @param A single adjacency matrix OR
#' @param A_by_time named list of adjacency matrices (names match unique(df[[time]]))
#' @param method "sbm" or "lsm"
#' @param directed whether adjacency is directed
#' @param add_homophily if TRUE, auto-add homophily terms to formula
#'
#' @param sbm_model SBM model: "bernoulli","poisson","gaussian" (used if method="sbm")
#' @param block_col name of SBM block column (used if method="sbm")
#' @param add_block_fe if TRUE, adds factor(block_col) to formula (SBM only)
#'
#' @param lsm_d embedding dimension (default 2) (LSM only)
#' @param lsm_cols prefix for LSM coordinate columns (LSM only)
#' @param add_lsm_linear if TRUE, adds the coordinate columns to formula (LSM only)
#'
#' @param fit_fun model fitting function (default stats::lm)
#' @param ... passed to fit_fun and/or to SBM/LSM detectors (as appropriate)
#'
#' @return list(model, data, formula, diagnostics, method_details)
#' @export
clean <- function(df, formula, id, time,
                  A = NULL,
                  A_by_time = NULL,
                  method = c("sbm","lsm"),
                  directed = TRUE,
                  add_homophily = TRUE,

                  # SBM options
                  sbm_model = c("bernoulli","poisson","gaussian"),
                  block_col = "clean_block",
                  add_block_fe = TRUE,

                  # LSM options
                  lsm_d = 2,
                  lsm_cols = "clean_lsm",
                  add_lsm_linear = TRUE,

                  # general
                  fit_fun = stats::lm,
                  seed = 1,
                  ...) {

  method <- match.arg(method)
  sbm_model <- match.arg(sbm_model)

  if (!is.data.frame(df)) stop("`df` must be a data.frame.")
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.")
  if (!all(c(id, time) %in% names(df))) stop("`id` and/or `time` not found in df.")
  if (!is.function(fit_fun)) stop("`fit_fun` must be a function.")

  # ---- decide network input mode ----
  if (is.null(A) && is.null(A_by_time)) {
    stop("Provide either `A` (single matrix) or `A_by_time` (named list of matrices).")
  }
  if (!is.null(A) && !is.null(A_by_time)) {
    stop("Provide only one of `A` or `A_by_time`, not both.")
  }

  # If user provides a single A, treat as one slice repeated for each time
  # (safe fallback; but we should message it clearly)
  if (!is.null(A)) {
    tt <- as.character(df[[time]])
    uu <- unique(tt)
    A_by_time <- stats::setNames(replicate(length(uu), A, simplify = FALSE), uu)
  }

  # ---- 1) augment df with homophily detection ----
  df2 <- df
  method_details <- list()

  if (method == "sbm") {

    # SBM: attach discrete block labels by time slice
    df2 <- clean_sbm_by_time(
      df = df2,
      id = id,
      time = time,
      A_by_time = A_by_time,
      model = sbm_model,
      directed = directed,
      seed = seed,
      ...
    )

    # clean_sbm_by_time() writes "clean_block" by default; we allow rename:
    # If user wants custom block_col, rename it.
    if (block_col != "clean_block") {
      if ("clean_block" %in% names(df2)) {
        df2[[block_col]] <- df2[["clean_block"]]
        df2[["clean_block"]] <- NULL
      } else if (!(block_col %in% names(df2))) {
        stop("Could not find SBM block column to rename; expected 'clean_block'.")
      }
    }

    method_details$block_col <- block_col
    method_details$model <- sbm_model

  } else {

    # LSM: should attach continuous coordinates by time slice.
    # EXPECTED interface of your function:
    # clean_lsm_by_time(df, id, time, A_by_time, d=2, directed=TRUE, out_prefix="clean_lsm", seed=1, ...)
    df2 <- clean_lsm_by_time(
      df = df2,
      id = id,
      time = time,
      A_by_time = A_by_time,
      d = lsm_d,
      directed = directed,
      out_prefix = lsm_cols,
      seed = seed,
      ...
    )

    method_details$d <- lsm_d
    method_details$lsm_cols <- paste0(lsm_cols, "_", seq_len(lsm_d))
  }

  # ---- 2) build augmented formula ----
  f2 <- formula
  ftxt <- paste(deparse(formula), collapse = "")

  if (add_homophily) {

    if (method == "sbm") {
      if (add_block_fe) {
        # add factor(block_col) if not already included
        if (!grepl(block_col, ftxt, fixed = TRUE)) {
          f2 <- stats::as.formula(paste0(ftxt, " + factor(", block_col, ")"))
        }
      }
    } else {
      # add linear coordinate controls (optional; users may prefer splines / interactions)
      if (add_lsm_linear) {
        coord_terms <- paste0(lsm_cols, "_", seq_len(lsm_d))
        missing_terms <- coord_terms[!vapply(coord_terms, function(z) grepl(z, ftxt, fixed = TRUE), logical(1))]
        if (length(missing_terms) > 0) {
          f2 <- stats::as.formula(paste0(ftxt, " + ", paste(missing_terms, collapse = " + ")))
        }
      }
    }
  }

  # ---- 3) fit model ----
  mod <- fit_fun(f2, data = df2, ...)

  # ---- 4) diagnostics ----
  diag <- list(
    n = nrow(df2),
    n_time = length(unique(as.character(df2[[time]]))),
    method = method,
    directed = directed
  )

  if (method == "sbm") {
    diag$block_col <- block_col
    diag$share_no_block <- mean(df2[[block_col]] == "no_block", na.rm = TRUE)
    diag$n_blocks <- length(unique(df2[[block_col]][df2[[block_col]] != "no_block"]))
  } else {
    coord_terms <- paste0(lsm_cols, "_", seq_len(lsm_d))
    diag$lsm_d <- lsm_d
    diag$missing_lsm <- any(!coord_terms %in% names(df2))
  }

  list(
    model = mod,
    data = df2,
    formula = f2,
    diagnostics = diag,
    method_details = method_details
  )
}
