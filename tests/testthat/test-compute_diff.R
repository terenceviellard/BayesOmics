#' @importFrom bayestestR overlap
#' @importFrom stats rnorm

test_that("calculate_group_overlaps computes the correct overlap matrix", {
  # Define the results for each group
  results <- list(
    group1 = list(muk = c(0), sigmak = matrix(c(1), nrow = 1)),
    group2 = list(muk = c(1), sigmak = matrix(c(1), nrow = 1)),
    group3 = list(muk = c(2), sigmak = matrix(c(1), nrow = 1))
  )

  # Define the number of samples
  n <- 10000

  # Calculate the overlap matrix
  overlap_matrix <- calculate_group_overlaps(results, n)

  # Check if the overlap matrix is symmetric
  expect_true(all(overlap_matrix == t(overlap_matrix)))

  # Check if the diagonal elements are all 1
  expect_true(all(diag(overlap_matrix) == 1))

  # Check if the overlap coefficients are within the expected range [0, 1]
  expect_true(all(overlap_matrix >= 0 & overlap_matrix <= 1))

  # Check specific overlap values if you have expectations
  # For example, the overlap between group1 and group2 should be less than 1
  expect_true(overlap_matrix["group1", "group2"] < 1)
  expect_true(overlap_matrix["group1", "group3"] < overlap_matrix["group1", "group2"])
})
