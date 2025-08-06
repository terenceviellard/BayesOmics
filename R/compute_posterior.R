#' Compute Posterior Means for Different Groups
#'
#' This function computes the posterior means for different groups within a dataset.
#' It expects the input data frame to contain a column named 'Group' that identifies
#' the groups, as well as columns 'ID', 'Output', and 'Input'.
#'
#' @param data A data frame containing the data to be analyzed. Must include columns 'Group', 'ID', 'Output', and 'Input'.
#' @param kern A kernel function or object used to compute pairwise kernels.
#' @param mu_0 Prior mean parameter.
#' @param lambda_0 Prior precision parameter.
#' @return A list of results for each group, containing the posterior means and kernel matrices.
#' @export
multi_posterior_mean <- function(data, kern, mu_0 = 1, lambda_0 = 1) {
  if (!"Group" %in% names(data)) {
    stop("The 'Group' column is not found in the input data.")
  }

  # Get unique groups and create group names
  groups <- unique(data$"Group")
  group_names <- paste0("Group", seq_along(groups))

  results <- list()

  for (i in seq_along(groups)) {
    group <- groups[i]
    group_name <- group_names[i]
    group_df <- subset(data, data$"Group" == group)
    all_peptides <- unique(group_df$ID)
    sum_outputs <- tapply(group_df$Output, group_df$ID, sum, na.rm = TRUE)
    lengths <- tapply(group_df$Output, group_df$ID, length)
    muk_vector <- (lambda_0 * mu_0 + sum_outputs) / (lengths + lambda_0)
    names(muk_vector) <- all_peptides

    if (nrow(group_df) > 0) {
      inputs <- as.matrix(group_df$Input)
      sigmak <- pairwise_kernel(kern, inputs, inputs) / (nrow(group_df) + lambda_0)
    } else {
      sigmak <- NA
    }

    results[[group_name]] <- list(muk = muk_vector, sigmak = sigmak)
  }

  return(results)
}





#' @title Sample from a Normal multivariate distribution
#'
#' @description
#' Sample
#'
#' @param results  A list containing parameters of the posterior distribution for each group.
#' Each group should have elements `muk` (mean) and `sigmak` (covariance matrix). (from multi_posterior_mean)
#' @param n A number indicating the number of samples
#'
#' @return A list of sampled values for each group.
#' @export
#'
sample_posterior <- function(results, n) {
  group_names <- names(results)
  samples_list <- list()

  for (group in group_names) {
    group_result <- results[[group]]

    # Check if muk and sigmak are present
    if (!all(c("muk", "sigmak") %in% names(group_result))) {
      stop("Each group in results must contain 'muk' and 'sigmak'.")
    }

    # Ensure sigmak is a matrix
    if (!is.matrix(group_result$sigmak)) {
      stop("sigmak must be a matrix.")
    }

    # Number of variables to sample
    num_vars <- length(group_result$muk)

    # Sample from the posterior distribution
    samples_matrix <- matrix(nrow = n, ncol = num_vars)
    for (i in 1:num_vars) {
      samples_matrix[, i] <- rnorm(n, mean = group_result$muk[i], sd = sqrt(group_result$sigmak[i, i]))
    }

    # Convert the matrix to a data frame
    samples_df <- as.data.frame(samples_matrix)
    names(samples_df) <- paste0("Var", 1:num_vars)

    # Store the data frame in the list
    samples_list[[group]] <- samples_df
  }

  return(samples_list)
}
