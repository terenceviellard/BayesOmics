# Generate a Synthetic Dataset for BayesOmics

Simulate a basic complete training dataset. Several flexible arguments
allow adjustment of the number of id, groups, and samples in each
experiment. The values of several parameters controlling the data
generation process can be modified.

## Usage

``` r
simu_db(
  nb_id = NULL,
  nb_group = 2,
  nb_sample = 5,
  nb_dim = 1,
  range_output = c(0, 50),
  range_input = c(0, 50),
  diff_group = 3,
  var_sample = 2,
  univariate = FALSE
)
```

## Arguments

- nb_id:

  An integer, indicating the number of id in the data. Ignored when
  `univariate = TRUE` (forced to `1`); leave at its default (`NULL`,
  resolved to `5` outside univariate mode and `1` inside it) unless
  overriding the multi-feature count.

- nb_group:

  An integer, indicating the number of groups/conditions.

- nb_sample:

  An integer, indicating the number of samples in the data for each id
  (i.e., the repetitions of the same experiment).

- nb_dim:

  An integer, indicating the number of Input dimensions per id. Defaults
  to `1` (a single scalar Input per id, emitted with `Input_ID = 1`);
  `nb_dim > 1` emits one row per (Group, ID, Sample, Input_ID), `Output`
  repeated identically across the `nb_dim` rows of one observation. Must
  stay `1` when `univariate = TRUE` (there is no Input axis to place in
  higher dimensions).

- range_output:

  A 2-sized vector, indicating the range of values for output from which
  to pick a mean value for each id

- range_input:

  A 2-sized vector, indicating the range of values for input from which
  to pick a mean value for each id (applied independently to every Input
  dimension). Ignored when `univariate = TRUE`.

- diff_group:

  A number, indicating the mean difference between consecutive groups.

- var_sample:

  A number, indicating the noise variance for each new sample of a id

- univariate:

  Logical. If `TRUE`, simulate a single feature per group (no
  feature/covariate axis) and return only `ID`, `Group`, `Sample`,
  `Output` – the format
  [`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
  requires for its univariate mode (see
  `dev/70_documentation/drafts/02_univariate.Rmd`). Requires `nb_id = 1`
  (or left at its default) and `nb_dim = 1` (or left at its default);
  `range_input` is unused. Defaults to `FALSE`.

## Value

A full dataset of synthetic data. With `univariate = FALSE` (default),
columns `ID`, `Group`, `Sample`, `Input_ID`, `Input`, `Output`. With
`univariate = TRUE`, columns `Group`, `Sample`, `Output` only (no
`ID`/`Input`) – ready to pass to
[`posterior_mean()`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
for its univariate mode.

## Examples

``` r
data <- simu_db()

# Univariate mode: a single feature per group, ready for posterior_mean()'s
# univariate mode (no ID/Input columns needed there):
data_uni <- simu_db(univariate = TRUE)
```
