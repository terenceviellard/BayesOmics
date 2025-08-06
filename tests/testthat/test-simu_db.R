test_that("simu_db generates correct dataset structure", {
  # Test with default parameters
  data <- simu_db()

  # Check if the output is a data frame
  expect_s3_class(data, "data.frame")

  # Check the number of rows
  expect_equal(nrow(data), 5 * 2 * 5) # nb_peptide * nb_group * nb_sample

  # Check column names
  expect_equal(colnames(data), c("ID", "Group", "Sample", "Input", "Output"))

  # Check if the Output column has the correct length
  expect_equal(nrow(data), length(data$Output))

  # Check if the Output values are within a reasonable range
  min_range <- 0
  max_range <- 50
  expect_true(all(data$Output >= min_range & data$Output <= max_range + 2 * 3 + 3 * 2))
})

test_that("simu_db respects custom parameters", {
  # Test with custom parameters
  custom_data <- simu_db(nb_id = 3, nb_group = 2, nb_sample = 4, range_output = c(10, 20), diff_group = 2, var_sample = 1)

  # Check if the output is a data frame
  expect_s3_class(custom_data, "data.frame")

  # Check the number of rows
  expect_equal(nrow(custom_data), 3 * 2 * 4) # nb_peptide * nb_group * nb_sample

  # Check if the Output values are within a reasonable range
  expect_true(all(custom_data$Output >= min(c(10, 20)) & custom_data$Output <= max(c(10, 20)) + 2 * 2 + 3 * 1))
})
