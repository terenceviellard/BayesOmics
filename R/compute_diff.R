#' @importFrom bayestestR overlap
#' @importFrom stats rnorm

#' @title Compute Overlapping Coefficient
#'
#' @description
#' Use the bayestestR package to compute the overlapping coefficient.
#'
#' @param results A list containing parameters of the posterior distribution for each group.
#' @param n A number indicating the number of samples generated during the overlapping coefficient computation process.
#'
#' @return A matrix with the coefficient for each group.
#' @export
#'
calculate_group_overlaps <- function(results, n = 1000) {
  group_names <- names(results)
  num_groups <- length(group_names)
  overlap_matrix <- matrix(NA, nrow = num_groups, ncol = num_groups)
  rownames(overlap_matrix) <- group_names
  colnames(overlap_matrix) <- group_names
  diag(overlap_matrix) <- 1
  for (i in 1:(num_groups - 1)) {
    for (j in (i + 1):num_groups) {
      group1 <- results[[i]]
      print(group1)
      group2 <- results[[j]]
      print(group2)
      # Check if muk and sigmak are present
      if (!all(c("muk", "sigmak") %in% names(group1)) || !all(c("muk", "sigmak") %in% names(group2))) {
        stop("Each group in results must contain 'muk' and 'sigmak'.")
      }

      sample1 <- rnorm(n, mean = group1$muk, sd = group1$sigmak)
      sample2 <- rnorm(n, mean = group2$muk, sd = group2$sigmak)
      overlap_coefficient <- overlap(sample1, sample2)
      overlap_matrix[group_names[i], group_names[j]] <- overlap_coefficient
      overlap_matrix[group_names[j], group_names[i]] <- overlap_coefficient
    }
  }
  return(overlap_matrix)
}
