#' @importFrom stats rnorm

# Test for ensure_abstract_kernel
test_that("ensure_abstract_kernel returns an AbstractKernel object", {
  # Test with a non-AbstractKernel input
  kernel <- ensure_abstract_kernel(1)
  expect_s4_class(kernel, "ConstantKernel")

  # Test with an AbstractKernel input
  abstract_kernel <- new("SEKernel")
  kernel <- ensure_abstract_kernel(abstract_kernel)
  expect_s4_class(kernel, "SEKernel")
})

# Test for get_hyperparameter_names
test_that("get_hyperparameter_names returns correct hyperparameter names", {
  kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5)
  hps <- get_hyperparameter_names(kernel)
  expect_equal(sort(hps), c("length_scale_se", "variance_se"))
})

# Test for get_hyperparameter_values
test_that("get_hyperparameter_values returns correct hyperparameter values", {
  kernel <- new("SEKernel", variance_se = 2.0, length_scale_se = 0.5)
  hps= get_hyperparameter_values(kernel)
  expect_equal(hps[["variance_se"]], 2)
  expect_equal(hps[["length_scale_se"]], 0.5)
})

# Test for set_hyperparameters
test_that("set_hyperparameters updates hyperparameters correctly", {
  kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5)
  updated_kernel <- set_hyperparameters(kernel, list(variance_se = 3, length_scale_se = 1))
  hps <- get_hyperparameter_values(updated_kernel)
  expect_equal(hps[["variance_se"]], 3)
  expect_equal(hps[["length_scale_se"]], 1)
})

# Test for chol_inv_jitter
test_that("chol_inv_jitter returns the correct inverse matrix", {
  mat <- matrix(c(2, 1, 1, 2), nrow = 2)
  expected_inv_mat <- matrix(c(0.6667,-0.3333, -0.3333, 0.6667), nrow = 2)
  inv_mat <- chol_inv_jitter(mat, pen_diag = 0.1)
  expect_equal(inv_mat, expected_inv_mat, tolerance = 1e-1)
})

# Test for dmnorm
test_that("dmnorm computes the correct density", {
  x <- matrix(c(1, 2, 3, 4), nrow = 2)
  mu <- c(0, 0)
  inv_Sigma <- diag(2)
  density <- dmnorm(x, mu, inv_Sigma, log = FALSE)
  expect_true(all(density > 0))
})


#test_that("optim_hp optimizes hyperparameters correctly with synthetic data", {
  # Define true hyperparameters
  #true_hp <- c(variance_se = 1, length_scale_se = 1)

  # Create a kernel with true hyperparameters
  #true_kernel <- new("SEKernel", variance_se = true_hp['variance_se'], length_scale_se = true_hp['length_scale_se'])

  # Generate synthetic input data
  #input <- seq(-5, 5, length.out = 100)
  #input <- matrix(input, ncol = 1)

  # Generate synthetic output data using the true kernel
  #cov_matrix <- pairwise_kernel(true_kernel, input, input)

  #output=stats::rnorm(1, mean = rep(0,100), sd = cov_matrix)

  # Create a data frame for the synthetic data
  #db <- data.frame(Input = input, Output = t(output))

  # Define the mean, kernel, and other parameters for optimization
  #mean <- 0
  #kern <- new("SEKernel", variance_se = 1.5, length_scale_se = 0.5) # Initial guess
  #post_cov <- 0.1
  #pen_diag <- 1e-4

  # Optimize the hyperparameters
  #optimized_hp <- optim_hp(c(1.1, 0.9), db, mean, kern, post_cov, pen_diag)

  # Define a tolerance level
  #tolerance <- 0.2

  # Print the optimized and true hyperparameters for debugging
  #print(optimized_hp)
  #print(true_hp)

  # Check if the optimized hyperparameters are close to the true values
  #expect_true(all(abs(optimized_hp - true_hp) < tolerance))
#})




# Test for ensure_abstract_kernel
test_that("ensure_abstract_kernel returns an AbstractKernel object", {
  # Test with a non-AbstractKernel input
  kernel <- ensure_abstract_kernel(1)
  expect_s4_class(kernel, "ConstantKernel")

  # Test with an AbstractKernel input
  abstract_kernel <- new("SEKernel")
  kernel <- ensure_abstract_kernel(abstract_kernel)
  expect_s4_class(kernel, "SEKernel")

  product_kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5) * new("SEKernel", variance_se = 1, length_scale_se = 1)
  ensured_product_kernel <- ensure_abstract_kernel(product_kernel)
  expect_s4_class(ensured_product_kernel, "ProductKernel")
})

# Test for get_hyperparameter_names
test_that("get_hyperparameter_names returns correct hyperparameter names", {
  kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5)
  hps <- get_hyperparameter_names(kernel)
  expect_equal(sort(hps), c("length_scale_se", "variance_se"))
})

test_that("get_hyperparameter_values returns correct hyperparameter values for product kernel", {
  # Create a product of SEKernel instances
  kernel <- new("SEKernel", variance_se = 2.0, length_scale_se = 0.5) * new("SEKernel", variance_se = 1.0, length_scale_se = 1.0)

  # Get hyperparameter values
  hps <- get_hyperparameter_values(kernel)

  expect_true(all(unlist(hps) == c(2.0, 0.5, 1.0, 1.0)))
})


test_that("set_hyperparameters updates hyperparameters correctly for product kernel", {
  # Create a product of SEKernel instances
  kernel <- new("SEKernel", variance_se = 2, length_scale_se = 0.5) * new("ConstantKernel", value_c = 1)
  updated_kernel <- set_hyperparameters(kernel, c(variance_se = 2.0, length_scale_se = 0.5, value_c = 2))

  # Get updated hyperparameter values
  hps <- get_hyperparameter_values(updated_kernel)
  expect_true(all(unlist(hps) == c(2, 0.5, 2)))
})


# Test for kernel operations
test_that("Kernel operations work correctly", {
  # Create instances of SEKernel
  K1 <- new("SEKernel", variance_se = 1, length_scale_se = 1)
  K2 <- new("SEKernel", variance_se = 2, length_scale_se = 0.5)
  K3 <- new("SEKernel", variance_se = 0.5, length_scale_se = 2)

  # Combine kernels using product and sum
  combined_kernel <- K1 * K2 + K3

  # Check if the combined kernel is of the correct class
  expect_s4_class(combined_kernel, "SumKernel")

  # Check if the kernels are combined correctly
  expect_true(length(combined_kernel@kernels) == 2)
  expect_s4_class(combined_kernel@kernels[[1]], "ProductKernel")
  expect_s4_class(combined_kernel@kernels[[2]], "SEKernel")
})



