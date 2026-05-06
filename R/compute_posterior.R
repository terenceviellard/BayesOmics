#' Compute Posterior Means for Different Groups
#'
#' This function computes the posterior means for different groups within a dataset.
#' It expects the input data frame to contain a column named 'Group' that identifies
#' the groups, as well as columns 'ID', 'Output', and 'Input'.
#'
#' The posterior distribution is given by:
#'
#' \eqn{p(\mathbf{\mu} \mid y_1, \dots, y_N, \Sigma_{\hat{\theta}}) = \mathcal{N}\left(\mathbf{\mu}; \ \dfrac{\lambda_0 \mu_0 + \sum_{n=1}^{N} y_n}{N + \lambda_0}, \dfrac{1}{N + \lambda_0} \Sigma_{\hat{\theta}}\right)}
#'
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom rlang .data
#' @param data A data frame containing the data to be analyzed. Must include columns 'Group', 'ID', 'Output', and 'Input'.
#' @param kern A kernel function or object used to compute pairwise kernels.
#' @param mu_0 Prior mean parameter.
#' @param lambda_0 Prior precision parameter.
#' @return A list of results for each group, containing the posterior means and kernel matrices.
#' @export
multi_posterior_mean <- function(data, kern, mu_0 = 1, lambda_0 = 1) {
  # === Initial checks ===
  required_cols <- c("Group", "ID", "Output", "Input")
  if (!all(required_cols %in% names(data))) {
    stop(paste0("The following columns are missing: ", paste(setdiff(required_cols, names(data)), collapse = ", ")))
  }
  if (!is.numeric(data$Input) || !is.numeric(data$Output) || !is.numeric(mu_0)) {
    stop("The 'Input' and 'Output' columns must be numeric, and mu_0 must be numeric.")
  }
  if (lambda_0 <= 0 || !is.numeric(lambda_0)) {
    stop("lambda_0 must be numeric and strictly positive.")
  }
  if (any(is.na(data$Group)) || any(is.na(data$ID))) {
    stop("The 'Group' and 'ID' columns must not contain NA values.")
  }
  if (!is(kern, "AbstractKernel")) {
    stop("The 'kern' argument must be a valid kernel object from the keRnel package.")
  }

  # === Convert Group to character if necessary ===
  if (!is.character(data$Group)) {
    data$Group <- as.character(data$Group)
  }

  # === Check for empty groups ===
  group_info <- data %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::summarise(
      n_id = dplyr::n_distinct(.data$ID),
      n_unique_inputs = dplyr::n_distinct(.data$Input),
      .groups = "drop"
    )

  empty_mask <- group_info$n_id == 0 | group_info$n_unique_inputs == 0
  if (any(empty_mask)) {
    empty_groups <- group_info$Group[empty_mask]
    stop(paste0("Empty groups: ", paste(empty_groups, collapse = ", ")))
  }

  # === Calculate df_mu ===
  df_mu <- data %>%
    dplyr::group_by(.data$Group, .data$ID) %>%
    dplyr::summarise(
      sum_outputs = sum(.data$Output, na.rm = TRUE),
      lengths = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(muk_vector = (lambda_0 * mu_0 + .data$sum_outputs) / (lengths + lambda_0))

  # === Calculate unique_inputs (list of sorted vectors) ===
  unique_inputs <- data %>%
    dplyr::group_by(.data$Group) %>%
    dplyr::summarise(unique_input = list(sort(unique(.data$Input))), .groups = "drop")

  # === Pre-compute subsets for group_data ===
  df_mu_by_group <- split(df_mu, df_mu$Group)
  unique_inputs_by_group <- stats::setNames(unique_inputs$unique_input, unique_inputs$Group)

  # === Build group_data ===
  groups <- unique(data$Group)
  group_data <- dplyr::tibble(
    group = groups,
    muk = purrr::map(groups, function(g) {
      subset <- df_mu_by_group[[g]]
      stats::setNames(subset$muk_vector, subset$ID)
    }),
    vec_name = purrr::map(groups, function(g) unique_inputs_by_group[[g]]),
    n_individus = purrr::map_int(groups, function(g) nrow(df_mu_by_group[[g]]))
  )

  # === Compute kernel matrices ===
  kernel_matrices <- list()  # Replaces new.env(hash = TRUE)

  results_list <- lapply(1:nrow(group_data), function(i) {
    vec <- group_data$vec_name[[i]]  # Already a vector
    vec_hash <- toString(sort(vec))  # Normalized key for caching
    if (!(vec_hash %in% names(kernel_matrices))) {
      kernel_matrices[[vec_hash]] <- keRnel::pairwise_kernel(
        kern,
        as.matrix(vec),
        as.matrix(vec)
      )
    }
    sigmak <- kernel_matrices[[vec_hash]] / (group_data$n_individus[i] + lambda_0)
    list(muk = group_data$muk[[i]], sigmak = sigmak)
  })
  names(results_list) <- group_data$group

  return(results_list)
}

#' @title Sample from a Normal multivariate distribution
#'
#' @description
#' Sample n elements from the posterior distribution.
#'
#' @importFrom mvtnorm rmvnorm
#' @param results  A list containing parameters of the posterior distribution for each group.
#' Each group should have elements `muk` (mean) and `sigmak` (covariance matrix). (from multi_posterior_mean)
#' @param n A number indicating the number of samples
#'
#' @return A list of sampled values for each group.
#' @export
#'
sample_posterior <- function(results, n) {
  if (!is.list(results)) {
    stop("'results' must be a list.")
  }
  if (length(results) == 0) {
    stop("'results' cannot be an empty list.")
  }
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }
  group_names <- names(results)
  samples_list <- list()

  for (group in group_names) {
    group_result <- results[[group]]

    # Check if muk and sigmak are present
    if (!all(c("muk", "sigmak") %in% names(group_result))) {
      stop("Each group in results must contain 'muk' and 'sigmak'.")
    }}

    # Ensure sigmak is a matrix
    if (!is.matrix(group_result$sigmak)) {
      stop("sigmak must be a matrix.")
    }

    samples_list <- lapply(results, function(group) {
      muk <- group$muk
      sigmak <- group$sigmak

      # Use rmvnorm from mvtnorm for multivariate sampling
      samples_matrix <- mvtnorm::rmvnorm(
        n = n,
        mean = muk,
        sigma = sigmak
      )

      # Set column names to the names of muk (IDs)
      colnames(samples_matrix) <- names(muk)

      # Convert to data frame
      as.data.frame(samples_matrix)
    })

    names(samples_list) <- names(results)

    return(samples_list)
}



