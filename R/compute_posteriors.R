multi_posterior_mean <- function(data, kern, mu_0, lambda_0) {

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
