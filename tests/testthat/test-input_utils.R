test_that("normalize_input_cols() is a no-op when Input_ID is already present", {
  db <- data.frame(ID = "a", Input_ID = 1, Input = 5)
  expect_identical(normalize_input_cols(db), db)
})

test_that("normalize_input_cols() injects Input_ID = 1 for a legacy scalar Input column", {
  db  <- data.frame(ID = c("a", "b"), Input = c(1, 2))
  out <- normalize_input_cols(db)
  expect_identical(out$Input_ID, c(1, 1))
  expect_identical(out$Input, db$Input)
})

test_that("normalize_input_cols() errors when neither Input_ID nor Input is present", {
  db <- data.frame(ID = "a", Output = 1)
  expect_error(normalize_input_cols(db), "must contain")
})

test_that("row_input_key() on a single-column matrix matches format_input_key()", {
  x   <- c(1.5, 2.25, -3)
  mat <- matrix(x, ncol = 1)
  expect_identical(row_input_key(mat), BayesOmics:::format_input_key(x))
})

test_that("row_input_key() coerces a plain vector to a 1-column matrix", {
  x <- c(1.5, 2.25, -3)
  expect_identical(row_input_key(x), row_input_key(matrix(x, ncol = 1)))
})

test_that("row_input_key() combines multiple columns into one key per row", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  keys <- row_input_key(mat)
  expect_length(keys, 2)
  expect_identical(
    keys,
    c(
      paste(BayesOmics:::format_input_key(mat[1, ]), collapse = "|"),
      paste(BayesOmics:::format_input_key(mat[2, ]), collapse = "|")
    )
  )
})

test_that("row_input_key() distinguishes rows differing in any dimension", {
  mat <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2, byrow = TRUE)
  expect_length(unique(row_input_key(mat)), 2)
})

test_that("input_matrix_by_id() pivots a D=1 long db into an N x 1 matrix", {
  db <- data.frame(
    ID = c("ID_1", "ID_2", "ID_3"),
    Input_ID = 1,
    Input = c(0.5, 3.1, -2)
  )
  pv <- input_matrix_by_id(db, key_cols = "ID")
  expect_equal(dim(pv$matrix), c(3, 1))
  expect_setequal(rownames(pv$matrix), db$ID)
  expect_equal(pv$matrix[db$ID, 1], setNames(db$Input, db$ID))
})

test_that("input_matrix_by_id() pivots a D=2 long db into an N x 2 matrix, column order by Input_ID", {
  db <- data.frame(
    ID       = rep(c("ID_1", "ID_2"), each = 2),
    Input_ID = rep(c("dim2", "dim1"), 2),
    Input    = c(3.1, 0.5, 4.2, 0.7)
  )
  pv <- input_matrix_by_id(db, key_cols = "ID")
  expect_equal(dim(pv$matrix), c(2, 2))
  # column order is sort(unique(Input_ID)) = c("dim1", "dim2")
  expect_equal(pv$matrix["ID_1", ], c(0.5, 3.1), ignore_attr = TRUE)
  expect_equal(pv$matrix["ID_2", ], c(0.7, 4.2), ignore_attr = TRUE)
})

test_that("input_matrix_by_id() is robust to a different row order per observation", {
  db1 <- data.frame(ID = "ID_1", Input_ID = c("a", "b"), Input = c(1, 2))
  db2 <- data.frame(ID = "ID_1", Input_ID = c("b", "a"), Input = c(2, 1))
  expect_equal(input_matrix_by_id(db1, "ID")$matrix, input_matrix_by_id(db2, "ID")$matrix)
})

test_that("input_matrix_by_id() errors on a duplicated Input_ID for the same observation", {
  db <- data.frame(ID = "ID_1", Input_ID = c(1, 1), Input = c(0.5, 0.6))
  expect_error(input_matrix_by_id(db, "ID"), "Duplicated")
})

test_that("input_matrix_by_id() duplicated-Input_ID error correctly names the offending ID and Input_ID (not NA)", {
  # Regression test: which(arr.ind = TRUE) labels its result columns after
  # table()'s own argument names, not literally "row"/"col" -- indexing by
  # those literal names silently produced NA/NA in the error message.
  db <- data.frame(ID = c("ID_1", "ID_1"), Input_ID = c(1, 1), Input = c(5, 6))
  expect_error(input_matrix_by_id(db, "ID"), "Duplicated 'Input_ID' value '1' for observation 'ID_1'")
})

test_that("input_matrix_by_id() errors on a ragged (inconsistent) Input_ID set across observations", {
  db <- data.frame(
    ID       = c("ID_1", "ID_1", "ID_2"),
    Input_ID = c(1, 2, 1),
    Input    = c(0.5, 3.1, 0.7)
  )
  expect_error(input_matrix_by_id(db, "ID"), "Ragged|Inconsistent")
})

test_that("input_matrix_by_id() supports multi-column key_cols", {
  db <- data.frame(
    Group    = rep(c("g1", "g2"), each = 2),
    ID       = rep(c("ID_1", "ID_2"), 2),
    Input_ID = 1,
    Input    = c(0.5, 3.1, 0.5, 3.1)
  )
  pv <- input_matrix_by_id(db, key_cols = c("Group", "ID"))
  expect_equal(nrow(pv$matrix), 4)
  expect_equal(nrow(pv$keys), 4)
})
