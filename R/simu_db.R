#' @noRd
validate_simu_db_common <- function(nb_id, nb_group, range_output, range_input, var_sample,
                                     var_sample_strictly_positive = FALSE) {
  if (!is.numeric(nb_id) || length(nb_id) != 1 || nb_id < 1 || nb_id != round(nb_id)) {
    stop("'nb_id' must be a single positive integer.")
  }
  if (!is.numeric(nb_group) || length(nb_group) != 1 || nb_group < 1 || nb_group != round(nb_group)) {
    stop("'nb_group' must be a single positive integer.")
  }
  if (!is.numeric(range_output) || length(range_output) != 2 || range_output[1] > range_output[2]) {
    stop("'range_output' must be a numeric vector of length 2 with range_output[1] <= range_output[2].")
  }
  if (!is.numeric(range_input) || length(range_input) != 2 || range_input[1] > range_input[2]) {
    stop("'range_input' must be a numeric vector of length 2 with range_input[1] <= range_input[2].")
  }
  if (var_sample_strictly_positive) {
    if (!is.numeric(var_sample) || length(var_sample) != 1 || var_sample <= 0) {
      stop("'var_sample' must be a single positive number.")
    }
  } else {
    if (!is.numeric(var_sample) || length(var_sample) != 1 || var_sample < 0) {
      stop("'var_sample' must be a single non-negative number.")
    }
  }
}

#' @noRd
validate_nb_dim <- function(nb_dim) {
  if (!is.numeric(nb_dim) || length(nb_dim) != 1 || nb_dim < 1 || nb_dim != round(nb_dim)) {
    stop("'nb_dim' must be a single positive integer.")
  }
}

#' @noRd
validate_simu_db_args <- function(nb_id, nb_group, nb_sample, range_output,
                                   range_input, diff_group, var_sample,
                                   var_sample_strictly_positive = FALSE) {
  validate_simu_db_common(nb_id, nb_group, range_output, range_input, var_sample,
                           var_sample_strictly_positive)
  if (!is.numeric(nb_sample) || length(nb_sample) != 1 || nb_sample < 1 || nb_sample != round(nb_sample)) {
    stop("'nb_sample' must be a single positive integer.")
  }
  if (!is.numeric(diff_group) || length(diff_group) != 1) {
    stop("'diff_group' must be a single number.")
  }
}

#' @importFrom stats rnorm runif

#' @title Generate a Synthetic Dataset for BayesOmics
#'
#' @description
#' Simulate a basic complete training dataset.
#' Several flexible arguments allow adjustment of the number of id, groups, and samples in each experiment.
#' The values of several parameters controlling the data generation process can be modified.
#'
#' @param nb_id An integer, indicating the number of id in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer, indicating the number of samples in the data for each id (i.e., the repetitions of the same experiment).
#' @param nb_dim An integer, indicating the number of Input dimensions per id. Defaults to `1` (a single scalar Input per id, emitted with `Input_ID = 1`); `nb_dim > 1` emits one row per (Group, ID, Sample, Input_ID), `Output` repeated identically across the `nb_dim` rows of one observation.
#' @param range_output A 2-sized vector, indicating the range of values for output from which to pick a mean value for each id
#' @param range_input A 2-sized vector, indicating the range of values for input from which to pick a mean value for each id (applied independently to every Input dimension)
#' @param diff_group A number, indicating the mean difference between consecutive groups.
#' @param var_sample A number, indicating the noise variance for each new sample of a id
#'
#' @return A full dataset of synthetic data, with columns `ID`, `Group`, `Sample`, `Input_ID`, `Input`, `Output`.
#' @export
#'
#' @examples
#' data <- simu_db()
simu_db <- function(
    nb_id = 5,
    nb_group = 2,
    nb_sample = 5,
    nb_dim = 1,
    range_output = c(0, 50),
    range_input = c(0, 50),
    diff_group = 3,
    var_sample = 2) {

  validate_simu_db_args(nb_id, nb_group, nb_sample, range_output, range_input,
                         diff_group, var_sample, var_sample_strictly_positive = FALSE)
  validate_nb_dim(nb_dim)

  base_output <- runif(nb_id, range_output[1], range_output[2])
  base_input  <- matrix(runif(nb_id * nb_dim, range_input[1], range_input[2]), nrow = nb_id, ncol = nb_dim)
  ids <- paste0("ID_", seq_len(nb_id))

  # One row per (ID, Group, Sample) -- the observation grain -- built once so
  # Output (and its noise draw) is computed a single time per observation,
  # then repeated identically across the nb_dim Input_ID rows below.
  obs <- data.frame(
    ID     = rep(ids, each = nb_group * nb_sample),
    Group  = rep(rep(seq_len(nb_group), each = nb_sample), nb_id),
    Sample = rep(seq_len(nb_sample), nb_group * nb_id),
    stringsAsFactors = FALSE
  )
  obs$Output <- rep(base_output, each = nb_group * nb_sample) +
               diff_group * (obs$Group - 1) +  # Group 1 = reference (effect 0)
               rnorm(nrow(obs), 0, var_sample)

  db <- do.call(rbind, lapply(seq_len(nb_dim), function(d) {
    cbind(obs, Input_ID = d, Input = rep(base_input[, d], each = nb_group * nb_sample))
  }))
  db <- db[, c("ID", "Group", "Sample", "Input_ID", "Input", "Output")]
  rownames(db) <- NULL

  return(db)
}


#' @importFrom stats rnorm runif
#' @importFrom methods is

#' @noRd
#'
#' @details Draws one column of `nb_id` Input positions (continuous, or
#'   integer grid/sample), independently of any other dimension -- used by
#'   `simu_db_kernel()` once per Input dimension, so the resulting `nb_id x
#'   nb_dim` matrix need not have globally-distinct D-dimensional tuples,
#'   only distinct values within each column.
simu_db_kernel_input_col <- function(nb_id, range_input, integer_input, input_grid) {
  if (!integer_input) {
    return(runif(nb_id, range_input[1], range_input[2]))
  }
  low  <- as.integer(ceiling(range_input[1]))
  high <- as.integer(floor(range_input[2]))
  if (input_grid) {
    # Regular grid: nb_id evenly spaced integer positions across [low, high]
    grid <- unique(round(seq(low, high, length.out = nb_id)))
    if (length(grid) < nb_id) {
      stop("Cannot generate ", nb_id, " distinct integers on a regular grid in [",
           low, ", ", high, "]; try a wider 'range_input'.")
    }
    return(as.numeric(sort(grid)))
  }
  # Random sampling without replacement from all integers in [low, high]
  pool <- seq(low, high)
  if (length(pool) < nb_id) {
    stop("'range_input' [", low, ", ", high, "] contains only ", length(pool),
         " distinct integer(s); cannot sample ", nb_id, " distinct values.")
  }
  as.numeric(sort(sample(pool, nb_id)))
}

#' @title Generate a Synthetic Dataset with Kernel-Structured Covariance
#'
#' @description
#' Simulate a complete training dataset, similar to \code{simu_db()}, but consistent
#' with the generative model underlying \code{\link{multi_posterior_mean}}: for each
#' group, the vector of Output values across the \code{nb_id} ids (indexed by their
#' Input position) is drawn jointly from a multivariate normal distribution whose
#' covariance is given by the kernel applied pairwise to the (shared, per-id) Input
#' values, plus independent measurement noise. Each replicate is an independent draw
#' of that same distribution:
#' \deqn{y_n \mid \mu_g \sim \mathcal{N}(\mu_g,\ \Sigma_\theta + \sigma^2 I), \quad n = 1,\dots,N_g}
#' where \eqn{\Sigma_\theta} is the kernel matrix over the (shared) Input values and
#' \eqn{\sigma^2} is \code{var_sample}. This makes \eqn{\Sigma_\theta} -- not just the
#' Input values -- shared across every group, exactly as assumed by
#' \code{\link{multi_posterior_mean}} and required by \code{\link{calculate_group_overlaps}}.
#'
#' By default (\code{mu_random = FALSE}) the per-group mean \eqn{\mu_g} is a
#' deterministic per-id baseline shifted by the corresponding element of \code{diff_group}.
#' Setting \code{mu_random = TRUE} instead draws \eqn{\mu_g} from the conjugate prior,
#' \eqn{\mu_g \sim \mathcal{N}(\mu_0,\ \Sigma_\theta/\lambda_0)} (plus the group shift)
#' -- useful for validating the posterior's frequentist coverage against the true \eqn{\mu_g}.
#'
#' @param nb_id An integer, indicating the number of features (ids) in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer or a vector of length \code{nb_group}. When a scalar,
#'   all groups receive the same number of replicates. When a vector, element \code{g}
#'   sets the number of replicates for group \code{g} independently, enabling unbalanced
#'   designs.
#' @param nb_dim An integer, indicating the number of Input dimensions per id.
#'   Defaults to \code{1} (a single scalar Input per id, emitted with
#'   \code{Input_ID = 1}). For \code{nb_dim > 1}, the kernel is evaluated
#'   directly on the \code{nb_id x nb_dim} matrix of positions (isotropic by
#'   default: a single shared \code{length_scale} across every dimension,
#'   unless \code{kernel} itself is an \code{ard_kernel()} or similar).
#'   \code{range_input} is applied independently to every dimension (each
#'   dimension's column of positions is drawn independently, not as
#'   globally-distinct D-dimensional tuples).
#' @param range_output A 2-element vector; its midpoint is used as the default baseline
#'   mean (\code{mu_0}) when \code{mu_0} is not supplied.
#' @param range_input A 2-element vector, indicating the range from which Input values
#'   are drawn for each feature and each dimension. Ignored for the integer dimension when \code{integer_input = TRUE}.
#' @param diff_group A numeric scalar or a vector of length \code{nb_group}. When a
#'   scalar, group \code{g} receives a shift of \code{diff_group * (g - 1)} (linear
#'   dose-response, group 1 = reference). When a vector, element \code{g} is the absolute
#'   offset applied to group \code{g}, allowing arbitrary non-linear group contrasts.
#' @param var_sample A positive number, the variance of the independent measurement noise
#'   added on top of the kernel-induced covariance.
#' @param kernel A kernel object (from the keRnel package). Defaults to
#'   \code{variance_kernel(variance = 1) * se_kernel(length_scale = 1)} when
#'   \code{NULL}.
#' @param mu_random If \code{TRUE}, draw each group's mean vector from the conjugate prior
#'   \eqn{\mathcal{N}(\mu_0, \Sigma_\theta/\lambda_0)} instead of a deterministic constant
#'   baseline. Defaults to \code{FALSE}.
#' @param mu_0 Baseline mean shared by every id (deterministic case) or prior mean
#'   (\code{mu_random = TRUE}); defaults to the midpoint of \code{range_output}.
#' @param lambda_0 Prior precision scaling (only used when \code{mu_random = TRUE});
#'   must be a single positive number. Defaults to \code{1}.
#' @param pen_diag Jitter added to the diagonal of the kernel-induced covariance matrix
#'   for numerical stability. Defaults to \code{1e-6}.
#' @param group_labels A character or numeric vector of length \code{nb_group} providing
#'   custom labels for the \code{Group} column. Defaults to \code{NULL}, which uses
#'   \code{"1"}, \code{"2"}, \ldots, \code{"nb_group"}.
#' @param integer_input Logical. If \code{TRUE}, Input values are distinct integers drawn
#'   from \code{[ceiling(range_input[1]), floor(range_input[2])]}. If \code{FALSE}
#'   (default), Input values are continuous uniform draws. See also \code{input_grid}.
#' @param input_grid Logical. Only used when \code{integer_input = TRUE}. If \code{FALSE}
#'   (default), integers are sampled uniformly without replacement (random positions). If
#'   \code{TRUE}, integers are placed on a regular grid of \code{nb_id} evenly spaced
#'   positions across the integer range (deterministic layout).
#'
#' @return A data frame of synthetic data with columns \code{ID}, \code{Group},
#'   \code{Sample}, \code{Input_ID}, \code{Input}, \code{Output}. Two attributes are attached:
#'   \describe{
#'     \item{\code{mu_true}}{Named list, one entry per group, each a named vector of true
#'       per-feature means.}
#'     \item{\code{base_input}}{An \code{nb_id x nb_dim} numeric matrix (\code{rownames = ID})
#'       giving the Input position(s) assigned to each feature.}
#'   }
#' @export
#'
#' @examples
#' # Default kernel (variance = 1, length_scale = 1):
#' data <- simu_db_kernel()
#'
#' # Custom kernel:
#' ker <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 3)
#' data <- simu_db_kernel(kernel = ker)
#'
#' # Unbalanced design: 3 replicates for group 1, 15 for group 2:
#' data <- simu_db_kernel(nb_sample = c(3, 15))
#'
#' # Arbitrary group contrasts (non-linear dose-response):
#' data <- simu_db_kernel(nb_group = 4, diff_group = c(0, 5, 5, 10))
#'
#' # Custom group labels:
#' data <- simu_db_kernel(nb_group = 2, group_labels = c("Control", "Treatment"))
#'
#' # Integer inputs on a regular grid:
#' data <- simu_db_kernel(nb_id = 10, integer_input = TRUE, input_grid = TRUE,
#'                        range_input = c(0, 100))
#'
#' # Multi-dimensional Input (D = 2), isotropic kernel:
#' data <- simu_db_kernel(nb_dim = 2)
simu_db_kernel <- function(
    nb_id = 5,
    nb_group = 2,
    nb_sample = 5,
    nb_dim = 1,
    range_output = c(0, 50),
    range_input = c(0, 50),
    diff_group = 3,
    var_sample = 2,
    kernel = NULL,
    mu_random = FALSE,
    mu_0 = NULL,
    lambda_0 = 1,
    pen_diag = 1e-6,
    group_labels = NULL,
    integer_input = FALSE,
    input_grid = FALSE) {

  # === Validation of common scalar parameters ===
  validate_simu_db_common(nb_id, nb_group, range_output, range_input, var_sample,
                           var_sample_strictly_positive = TRUE)
  validate_nb_dim(nb_dim)

  # nb_sample: scalar (same for all groups) or vector (one value per group)
  if (!is.numeric(nb_sample) || any(!is.finite(nb_sample)) ||
      any(nb_sample < 1) || any(nb_sample != round(nb_sample))) {
    stop("'nb_sample' must be a positive integer or a vector of positive integers.")
  }
  if (length(nb_sample) == 1) {
    nb_sample <- rep(as.integer(nb_sample), nb_group)
  } else if (length(nb_sample) == nb_group) {
    nb_sample <- as.integer(nb_sample)
  } else {
    stop("'nb_sample' must be a scalar or a vector of length 'nb_group'.")
  }

  # diff_group: scalar -> linear dose-response; vector -> explicit per-group offsets
  if (!is.numeric(diff_group) || any(!is.finite(diff_group))) {
    stop("'diff_group' must be a numeric scalar or a vector of length 'nb_group'.")
  }
  if (length(diff_group) == 1) {
    diff_group <- diff_group * (seq_len(nb_group) - 1)
  } else if (length(diff_group) == nb_group) {
    diff_group <- as.numeric(diff_group)
  } else {
    stop("'diff_group' must be a scalar or a vector of length 'nb_group'.")
  }

  # kernel: NULL -> default variance_kernel(1) * se_kernel(1)
  if (is.null(kernel)) {
    kernel <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
  } else if (!inherits(kernel, "kernel")) {
    stop("'kernel' must be a valid kernel object from the keRnel package (inheriting from 'kernel').")
  }

  if (!is.numeric(pen_diag) || length(pen_diag) != 1 || pen_diag < 0) {
    stop("'pen_diag' must be a single non-negative number.")
  }
  if (!is.logical(mu_random) || length(mu_random) != 1 || is.na(mu_random)) {
    stop("'mu_random' must be a single TRUE/FALSE value.")
  }
  if (mu_random && (!is.numeric(lambda_0) || length(lambda_0) != 1 || lambda_0 <= 0)) {
    stop("'lambda_0' must be a single positive number when mu_random = TRUE.")
  }
  if (!is.logical(integer_input) || length(integer_input) != 1 || is.na(integer_input)) {
    stop("'integer_input' must be a single TRUE/FALSE value.")
  }
  if (!is.logical(input_grid) || length(input_grid) != 1 || is.na(input_grid)) {
    stop("'input_grid' must be a single TRUE/FALSE value.")
  }

  # group_labels: NULL -> "1", "2", ...; otherwise must be unique and of length nb_group
  if (is.null(group_labels)) {
    labels <- as.character(seq_len(nb_group))
  } else {
    if (length(group_labels) != nb_group) {
      stop("'group_labels' must have length 'nb_group'.")
    }
    if (any(duplicated(as.character(group_labels)))) {
      stop("'group_labels' must contain unique values.")
    }
    labels <- as.character(group_labels)
  }

  # === Input generation ===
  # One column per Input dimension, shared across every group, so the kernel matrix
  # is the same for every group -- required by multi_posterior_mean() and
  # calculate_group_overlaps(). Each dimension is drawn independently (its own
  # distinct-integer draw/grid within the shared range_input, when
  # integer_input = TRUE) rather than requiring globally-distinct D-tuples.
  # matrix(vapply(...), nrow = nb_id) rather than relying on vapply's own
  # simplification: with nb_id = 1, FUN.VALUE has length 1, and vapply then
  # returns a bare vector (never a matrix) regardless of nb_dim.
  base_input <- matrix(
    vapply(
      seq_len(nb_dim),
      function(d) simu_db_kernel_input_col(nb_id, range_input, integer_input, input_grid),
      numeric(nb_id)
    ),
    nrow = nb_id, ncol = nb_dim
  )
  ids <- paste0("ID_", seq_len(nb_id))
  rownames(base_input) <- ids

  # === Kernel matrix ===
  Sigma_theta <- keRnel::evaluate(kernel, base_input, base_input)

  # Per-replicate observation covariance: kernel structure + independent measurement noise
  Sigma_obs <- chol_inv_jitter_diag(Sigma_theta + var_sample * diag(nb_id), pen_diag)
  L_obs     <- t(base::chol(Sigma_obs))

  mu_0_val <- if (is.null(mu_0)) mean(range_output) else mu_0
  if (mu_random) {
    L_prior <- t(base::chol(chol_inv_jitter_diag(Sigma_theta / lambda_0, pen_diag)))
  }

  # === Simulate one data frame per group ===
  groups <- lapply(seq_len(nb_group), function(g) {
    mu_g <- if (mu_random) {
      as.vector(mu_0_val + L_prior %*% rnorm(nb_id))
    } else {
      rep(mu_0_val, nb_id)
    } + diff_group[g]

    ns <- nb_sample[g]
    df <- do.call(rbind, lapply(seq_len(ns), function(s) {
      y   <- mu_g + as.vector(L_obs %*% rnorm(nb_id))
      obs <- data.frame(ID = ids, Group = labels[g], Sample = s, Output = y, stringsAsFactors = FALSE)
      do.call(rbind, lapply(seq_len(nb_dim), function(d) {
        cbind(obs, Input_ID = d, Input = base_input[, d])
      }))
    }))
    list(mu_true = stats::setNames(mu_g, ids), data = df)
  })

  result  <- do.call(rbind, lapply(groups, `[[`, "data"))
  result  <- result[, c("ID", "Group", "Sample", "Input_ID", "Input", "Output")]
  rownames(result) <- NULL
  mu_true <- lapply(groups, `[[`, "mu_true")
  attr(result, "mu_true")    <- stats::setNames(mu_true, labels)
  attr(result, "base_input") <- base_input
  result
}

#' @noRd
chol_inv_jitter_diag <- function(mat, pen_diag, max_tries = 20, warn_ratio = 100) {
  jitter_until_pd(mat, pen_diag, function(m) { chol(m); m },
                  max_tries, warn_ratio, label = "chol_inv_jitter_diag")
}
