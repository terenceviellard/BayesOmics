test_that("multi_posterior_mean computes posterior means correctly", {
  # Create a sample dataset
  data <- data.frame(
    Group = c(rep("Group1", 5), rep("Group2", 5)),
    ID = c(rep("ID1", 3), rep("ID2", 2), rep("ID3", 3), rep("ID4", 2)),
    Output = rnorm(10),
    Input = rnorm(10)
  )

  # Define a simple kernel function for testing
  kernel = new("SEKernel")
  # Parameters for the prior distribution
  mu_0 <- 0
  lambda_0 <- 1

  # Compute posterior means
  results <- multi_posterior_mean(data, kernel, mu_0, lambda_0)

  # Check if the results list contains entries for each group
  expect_equal(names(results), unique(data$Group))

  # Check if each group result contains 'muk' and 'sigmak'
  for (group in names(results)) {
    expect_true(all(c("muk", "sigmak") %in% names(results[[group]])))
    expect_equal(names(results[[group]]$muk), unique(data[data$Group == group, "ID"]))
  }
})

