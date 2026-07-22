#' @noRd
#'
#' @details Normalizes a `db` data frame to the long multi-dimensional Input
#'   format (`Input_ID`, `Input`) expected by every internal function past
#'   this point. A legacy `db` with only a scalar `Input` column (no
#'   `Input_ID`) gets `Input_ID = 1` injected -- this is the single
#'   normalization point that keeps the rest of the package (compute_posterior.R,
#'   optim_kernel.R, simu_db.R) working on one code path regardless of which
#'   format the caller supplied. `db` with `Input_ID` already present is
#'   returned unchanged (assumed already long-format, multi-dimensional or not).
normalize_input_cols <- function(db) {
  if ("Input_ID" %in% names(db)) {
    return(db)
  }
  if ("Input" %in% names(db)) {
    db$Input_ID <- 1
    return(db)
  }
  stop(
    "'db' must contain either 'Input_ID' and 'Input' columns (long, ",
    "possibly multi-dimensional format) or a single 'Input' column ",
    "(legacy scalar format)."
  )
}

#' @noRd
#'
#' @details One deterministic string key per row of an N x D input matrix,
#'   built from format_input_key()'s existing round-trip-safe "%.17g"
#'   formatting (compute_posterior.R). For D = 1 this produces output
#'   identical to `format_input_key(x)` applied directly to the column
#'   vector, so every existing scalar-Input caller sees no behavior change.
#'   This is the single primitive used to key/dedupe/align input points by
#'   value, both for the kernel-matrix cache in compute_posterior.R and the
#'   replicate-alignment logic in optim_kernel.R.
row_input_key <- function(mat) {
  if (!is.matrix(mat)) {
    mat <- matrix(mat, ncol = 1)
  }
  apply(mat, 1, function(row) paste(format_input_key(row), collapse = "|"))
}

#' @noRd
#'
#' @details Pivots a `db` already normalized by `normalize_input_cols()`
#'   (long format, `Input_ID`/`Input` columns) into an N x D numeric matrix:
#'   one row per unique combination of `key_cols` (e.g. `"ID"`), one column
#'   per unique value of `Input_ID`. Column order is fixed once, by
#'   `sort(unique(db$Input_ID))`, so the same physical dimension always lands
#'   in the same matrix column regardless of row order in `db` -- this is
#'   what keeps dimension identity aligned when the result is later compared
#'   or combined with another such matrix (e.g. across groups in
#'   `multi_posterior_mean()`, or across replicates in `optim_kernel.R`).
#'
#'   Errors if any `key_cols` combination has more than one row for the same
#'   `Input_ID` (a duplicated dimension for one observation), or if the set
#'   of `Input_ID` values is not identical across every `key_cols`
#'   combination (a "ragged" multi-dimensional design, unsupported).
#'
#' @param db A data frame already passed through `normalize_input_cols()`.
#' @param key_cols Character vector of column names that jointly identify one
#'   observation (e.g. `"ID"`, or `c("Group", "Sample", "ID")`).
#' @return A list with `matrix` (an N x D numeric matrix, `rownames` set to a
#'   `"_"`-joined key of `key_cols` when `length(key_cols) > 1`, or the bare
#'   `key_cols` values when `length(key_cols) == 1`) and `keys` (a data frame,
#'   `key_cols` columns only, same row order as `matrix`).
input_matrix_by_id <- function(db, key_cols) {
  if (!is.atomic(db$Input_ID) || anyNA(db$Input_ID)) {
    stop("'Input_ID' must be an atomic column with no missing values.")
  }
  dims <- sort(unique(db$Input_ID))
  d <- length(dims)

  row_key <- function(x) do.call(paste, c(x[key_cols], list(sep = "_")))

  keys_df <- unique(db[key_cols])
  key_strings <- row_key(keys_df)
  ord <- order(key_strings)
  keys_df <- keys_df[ord, , drop = FALSE]
  key_strings <- key_strings[ord]
  n <- length(key_strings)

  db_row_id <- match(row_key(db), key_strings)

  occ_counts <- table(db_row_id, db$Input_ID)
  if (any(occ_counts > 1)) {
    # Positional indexing (not by "row"/"col" name): which(arr.ind = TRUE)
    # labels its result columns after table()'s own argument names
    # (`db_row_id`, `db$Input_ID` here), not literally "row"/"col".
    bad <- which(occ_counts > 1, arr.ind = TRUE)[1, ]
    stop(
      "Duplicated 'Input_ID' value '", colnames(occ_counts)[bad[2]],
      "' for observation '", key_strings[as.integer(rownames(occ_counts)[bad[1]])],
      "' (key columns: ", paste(key_cols, collapse = ", "), ")."
    )
  }

  observed_dims <- split(db$Input_ID, db_row_id)
  ragged <- vapply(observed_dims, function(x) !setequal(x, dims), logical(1))
  if (any(ragged)) {
    bad_row_id <- as.integer(names(observed_dims)[ragged][1])
    stop(
      "Inconsistent 'Input_ID' set across observations: observation '",
      key_strings[bad_row_id], "' has Input_ID = {",
      paste(sort(unique(observed_dims[[as.character(bad_row_id)]])), collapse = ", "),
      "}, expected {", paste(dims, collapse = ", "), "}. Ragged multi-dimensional ",
      "designs (a different set of dimensions per observation) are not supported."
    )
  }

  mat <- matrix(NA_real_, nrow = n, ncol = d)
  mat[cbind(db_row_id, match(db$Input_ID, dims))] <- db$Input
  rownames(mat) <- key_strings

  list(matrix = mat, keys = keys_df)
}
