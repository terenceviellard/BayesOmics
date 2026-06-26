#' @noRd
validate_simu_db_args <- function(nb_id, nb_group, nb_sample, range_output,
                                   range_input, diff_group, var_sample,
                                   var_sample_strictly_positive = FALSE) {
  if (!is.numeric(nb_id) || length(nb_id) != 1 || nb_id < 1 || nb_id != round(nb_id)) {
    stop("'nb_id' must be a single positive integer.")
  }
  if (!is.numeric(nb_group) || length(nb_group) != 1 || nb_group < 1 || nb_group != round(nb_group)) {
    stop("'nb_group' must be a single positive integer.")
  }
  if (!is.numeric(nb_sample) || length(nb_sample) != 1 || nb_sample < 1 || nb_sample != round(nb_sample)) {
    stop("'nb_sample' must be a single positive integer.")
  }
  if (!is.numeric(range_output) || length(range_output) != 2 || range_output[1] > range_output[2]) {
    stop("'range_output' must be a numeric vector of length 2 with range_output[1] <= range_output[2].")
  }
  if (!is.numeric(range_input) || length(range_input) != 2 || range_input[1] > range_input[2]) {
    stop("'range_input' must be a numeric vector of length 2 with range_input[1] <= range_input[2].")
  }
  if (!is.numeric(diff_group) || length(diff_group) != 1) {
    stop("'diff_group' must be a single number.")
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

#' @importFrom stats rnorm runif

#' @title Generate a Synthetic Dataset Tailored for ProteoBayes
#'
#' @description
#' Simulate a basic complete training dataset.
#' Several flexible arguments allow adjustment of the number of id, groups, and samples in each experiment.
#' The values of several parameters controlling the data generation process can be modified.
#'
#' @param nb_id An integer, indicating the number of id in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer, indicating the number of samples in the data for each id (i.e., the repetitions of the same experiment).
#' @param range_output A 2-sized vector, indicating the range of values for output from which to pick a mean value for each id
#' @param range_input A 2-sized vector, indicating the range of values for input from which to pick a mean value for each id
#' @param diff_group A number, indicating the mean difference between consecutive groups.
#' @param var_sample A number, indicating the noise variance for each new sample of a id
#'
#' @return A full dataset of synthetic data.
#' @export
#'
#' @examples
#' data <- simu_db()
simu_db <- function(
    nb_id = 5,
    nb_group = 2,
    nb_sample = 5,
    range_output = c(0, 50),
    range_input = c(0, 50),
    diff_group = 3,
    var_sample = 2) {

  validate_simu_db_args(nb_id, nb_group, nb_sample, range_output, range_input,
                         diff_group, var_sample, var_sample_strictly_positive = FALSE)

  base_output <- runif(nb_id, range_output[1], range_output[2])
  base_input  <- runif(nb_id, range_input[1], range_input[2])

  db <- data.frame(
    ID     = rep(paste0("ID_", seq_len(nb_id)), each = nb_group * nb_sample),
    Group  = rep(rep(seq_len(nb_group), each = nb_sample), nb_id),
    Sample = rep(seq_len(nb_sample), nb_group * nb_id),
    Input  = rep(base_input, each = nb_group * nb_sample),
    stringsAsFactors = FALSE
  )

  db$Output <- rep(base_output, each = nb_group * nb_sample) +
               diff_group * (db$Group - 1) +  # Group 1 = reference (effect 0)
               rnorm(nrow(db), 0, var_sample)

  return(db)
}


#' @importFrom stats rnorm runif
#' @importFrom methods is

#' @title Generate a Synthetic Dataset with Kernel-Structured Covariance
#'
#' @description
#' Simulate a complete training dataset, similar to \code{simu_db()}, but consistent
#' with the generative model underlying \code{\link{multi_posterior_mean}}: for each
#' group, the vector of Output values across the \code{nb_id} ids (indexed by their
#' Input position) is drawn jointly from a multivariate normal distribution whose
#' covariance is given by the kernel applied pairwise to the (shared, per-id) Input
#' values, plus independent measurement noise. Each of the \code{nb_sample} replicates
#' of a group is an independent draw of that same distribution:
#' \deqn{y_n \mid \mu_g \sim \mathcal{N}(\mu_g,\ \Sigma_\theta + \sigma^2 I), \quad n = 1,\dots,\code{nb_sample}}
#' where \eqn{\Sigma_\theta} is the kernel matrix over the (shared) Input values and
#' \eqn{\sigma^2} is \code{var_sample}. This makes \eqn{\Sigma_\theta} -- not just the
#' Input values -- shared across every group, exactly as assumed by
#' \code{\link{multi_posterior_mean}} and required by \code{\link{calculate_group_overlaps}}.
#'
#' By default (\code{mu_random = FALSE}) the per-group mean \eqn{\mu_g} is a
#' deterministic per-id baseline shifted by \code{diff_group} per group level. Setting
#' \code{mu_random = TRUE} instead draws \eqn{\mu_g} from the conjugate prior assumed by
#' \code{\link{multi_posterior_mean}}, \eqn{\mu_g \sim \mathcal{N}(\mu_0,\ \Sigma_\theta/\lambda_0)}
#' (shifted by \code{diff_group} per group level) -- useful for validating the posterior's
#' frequentist coverage against the true generative \eqn{\mu_g}.
#'
#' @param nb_id An integer, indicating the number of id in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer, indicating the number of independent replicates per group.
#' @param range_output A 2-sized vector; its midpoint is used as the default baseline mean (\code{mu_0}) when \code{mu_0} is not supplied.
#' @param range_input A 2-sized vector, indicating the range of values for input from which to pick a position for each id.
#' @param diff_group A number, indicating the mean difference between consecutive groups.
#' @param var_sample A positive number, the variance of the independent measurement noise added on top of the kernel-induced covariance.
#' @param kernel A kernel object inheriting from \code{AbstractKernel} (keRnel package), used through \code{keRnel::pairwise_kernel()} to build the covariance matrix shared by every group.
#' @param mu_random If \code{TRUE}, draw each group's mean vector from the conjugate prior \eqn{\mathcal{N}(\mu_0, \Sigma_\theta/\lambda_0)} instead of using a deterministic constant baseline. Defaults to \code{FALSE}.
#' @param mu_0 Baseline mean shared by every id (deterministic case) or prior mean (\code{mu_random = TRUE}); defaults to the midpoint of \code{range_output}.
#' @param lambda_0 Prior precision scaling (only used when \code{mu_random = TRUE}); must be a single positive number. Defaults to \code{1}.
#' @param pen_diag Jitter added to the diagonal of the kernel-induced covariance matrix for numerical stability. Defaults to \code{1e-6}.
#'
#' @return A full dataset of synthetic data.
#' @export
#'
#' @examples
#' ker <- methods::new("SEKernel")
#' data <- simu_db_kernel(kernel = ker)
simu_db_kernel <- function(
    nb_id = 5,
    nb_group = 2,
    nb_sample = 5,
    range_output = c(0, 50),
    range_input = c(0, 50),
    diff_group = 3,
    var_sample = 2,
    kernel,
    mu_random = FALSE,
    mu_0 = NULL,
    lambda_0 = 1,
    pen_diag = 1e-6) {

  # === Input checks ===
  validate_simu_db_args(nb_id, nb_group, nb_sample, range_output, range_input,
                         diff_group, var_sample, var_sample_strictly_positive = TRUE)
  if (missing(kernel) || !methods::is(kernel, "AbstractKernel")) {
    stop("'kernel' must be a valid kernel object from the keRnel package (inheriting from 'AbstractKernel').")
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

  # One input value per id, shared across every group, so that the kernel matrix
  # (Sigma_theta below) is the SAME for every group -- mirroring multi_posterior_mean()'s
  # assumption and letting calculate_group_overlaps() compare these groups directly.
  base_input  <- runif(nb_id, range_input[1], range_input[2])
  input_mat   <- as.matrix(base_input)
  Sigma_theta <- keRnel::pairwise_kernel(kernel, input_mat, input_mat)

  # Per-replicate observation covariance: the kernel structure across ids, plus
  # independent measurement noise (var_sample), i.e. y_n | mu_g ~ N(mu_g, Sigma_theta + var_sample*I).
  Sigma_obs <- chol_inv_jitter_diag(Sigma_theta + var_sample * diag(nb_id), pen_diag)
  L_obs     <- t(base::chol(Sigma_obs))

  mu_0_val <- if (is.null(mu_0)) mean(range_output) else mu_0
  if (mu_random) {
    L_prior <- t(base::chol(chol_inv_jitter_diag(Sigma_theta / lambda_0, pen_diag)))
  }

  groups <- lapply(seq_len(nb_group), function(g) {
    mu_g <- if (mu_random) {
      as.vector(mu_0_val + L_prior %*% rnorm(nb_id))
    } else {
      rep(mu_0_val, nb_id)
    } + diff_group * (g - 1)  # Group 1 = reference (effect 0)

    df <- do.call(rbind, lapply(seq_len(nb_sample), function(s) {
      data.frame(
        ID     = paste0("ID_", seq_len(nb_id)),
        Group  = g,
        Sample = s,
        Input  = base_input,
        Output = mu_g + as.vector(L_obs %*% rnorm(nb_id)),
        stringsAsFactors = FALSE
      )
    }))
    list(mu_true = stats::setNames(mu_g, paste0("ID_", seq_len(nb_id))), data = df)
  })

  result <- do.call(rbind, lapply(groups, `[[`, "data"))

  mu_true <- lapply(groups, `[[`, "mu_true")
  attr(result, "mu_true") <- stats::setNames(mu_true, as.character(seq_len(nb_group)))
  attr(result, "base_input") <- stats::setNames(base_input, paste0("ID_", seq_len(nb_id)))
  result
}

#' @noRd
chol_inv_jitter_diag <- function(mat, pen_diag, max_tries = 20, warn_ratio = 100) {
  jitter_until_pd(mat, pen_diag, function(m) { chol(m); m },
                  max_tries, warn_ratio, label = "chol_inv_jitter_diag")
}
