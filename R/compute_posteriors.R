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
multi_posterior_mean <- function(data, kern, mu_0, lambda_0) {

  if (!"Group" %in% names(data)) {
    stop("The 'Group' column is not found in the input data.")
  }

  groups <- unique(data$Group)
  results <- list()

  for (group in groups) {
    group_df <- subset(data, Group == group)
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
    results[[group]] <- list(muk = muk_vector, sigmak = sigmak)
  }
  return(results)
}
