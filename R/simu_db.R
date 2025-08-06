#' @importFrom stats rnorm runif

#' @title Generate a Synthetic Dataset Tailored for ProteoBayes
#'
#' @description
#' Simulate a complete training dataset, which may be representative of various applications.
#' Several flexible arguments allow adjustment of the number of peptides, groups, and samples in each experiment.
#' The values of several parameters controlling the data generation process can be modified.
#'
#' @param nb_peptide An integer, indicating the number of peptides in the data.
#' @param nb_group An integer, indicating the number of groups/conditions.
#' @param nb_sample An integer, indicating the number of samples in the data for each peptide (i.e., the repetitions of the same experiment).
#' @param range_peptide A 2-sized vector, indicating the range of values from which to pick a mean value for each peptide.
#' @param diff_group A number, indicating the mean difference between consecutive groups.
#' @param var_sample A number, indicating the noise variance for each new sample of a peptide.
#'
#' @return A full dataset of synthetic data.
#' @export
#'
#' @examples
#' data <- simu_db(nb_peptide = 5, nb_group = 2, nb_sample = 3)
simu_db <- function(
    nb_peptide = 5,
    nb_group = 2,
    nb_sample = 5,
    range_peptide = c(0, 50),
    diff_group = 3,
    var_sample = 2
) {
  db <- data.frame(
    Peptide = rep(paste0('Peptide_', 1:nb_peptide), each = nb_group * nb_sample),
    Group = rep(rep(1:nb_group, each = nb_sample), nb_peptide),
    Sample = rep(1:nb_sample, nb_group * nb_peptide),
    stringsAsFactors = FALSE
  )

  db$Output <- unlist(by(db, db$Peptide, function(x) {
    runif(1, range_peptide[1], range_peptide[2])
  }))

  db$Output <- unlist(by(db, db$Group, function(x) {
    x$Output + diff_group * x$Group[1]
  }))

  db$Output <- db$Output + rnorm(nrow(db), 0, var_sample)

  return(db)
}
