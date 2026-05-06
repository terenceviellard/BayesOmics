#' @importFrom stats rnorm runif

#' @title Generate a Synthetic Dataset Tailored for ProteoBayes
#'
#' @description
#' Simulate a complete training dataset, which may be representative of various applications.
#' Several flexible arguments allow adjustment of the number of id, groups, and samples in each experiment.
#' The values of several parameters controlling the data generation process can be modified.
#'
#' @param nb_id An integer, indicating the number of id in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer, indicating the number of samples in the data for each id (i.e., the repetitions of the same experiment).
#' @param range_output A 2-sized vector, indicating the range of values for output from which to pick a mean value for each id
#' @param range_input A 2-sized vector, indicating the range of values for input from which to pick a mean value for each id
#' @param diff_group A number, indicating the mean difference between consecutive groups.
#' @param var_sample A number, indicating the noise variance for each new sample of a id
#'
#' @return A full dataset of synthetic data.
#' @export
#'
#' @examples
#' data <- simu_db()
simu_db <- function(
    nb_id = 5,
    nb_group = 2,
    nb_sample = 5,
    range_output = c(0, 50),
    range_input = c(0, 50),
    diff_group = 3,
    var_sample = 2) {
  db <- data.frame(
    ID = rep(paste0("ID_", 1:nb_id), each = nb_group * nb_sample),
    Group = rep(rep(1:nb_group, each = nb_sample), nb_id),
    Sample = rep(1:nb_sample, nb_group * nb_id),
    Input = runif(nb_group * nb_id, range_input[1], range_input[2]),
    stringsAsFactors = FALSE
  )
  # Generate Output column
  db$Output <- unlist(by(db, db$ID, function(x) {
    runif(1, range_output[1], range_output[2])
  }))
  db$Output <- unlist(by(db, db$Group, function(x) {
    x$Output + diff_group * x$Group[1]
  }))
  db$Output <- db$Output + rnorm(nrow(db), 0, var_sample)

  return(db)
}


