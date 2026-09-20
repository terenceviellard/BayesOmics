# Generate a Synthetic Dataset for BayesOmics

Simulate a basic complete training dataset. Several flexible arguments
allow adjustment of the number of id, groups, and samples in each
experiment. The values of several parameters controlling the data
generation process can be modified.

## Usage

``` r
simu_db(
  nb_id = 5,
  nb_group = 2,
  nb_sample = 5,
  nb_dim = 1,
  range_output = c(0, 50),
  range_input = c(0, 50),
  diff_group = 3,
  var_sample = 2
)
```

## Arguments

- nb_id:

  An integer, indicating the number of id in the data.

- nb_group:

  An integer, indicating the number of groups/conditions.

- nb_sample:

  An integer, indicating the number of samples in the data for each id
  (i.e., the repetitions of the same experiment).

- nb_dim:

  An integer, indicating the number of Input dimensions per id. Defaults
  to `1` (a single scalar Input per id, emitted with `Input_ID = 1`);
  `nb_dim > 1` emits one row per (Group, ID, Sample, Input_ID), `Output`
  repeated identically across the `nb_dim` rows of one observation.

- range_output:

  A 2-sized vector, indicating the range of values for output from which
  to pick a mean value for each id

- range_input:

  A 2-sized vector, indicating the range of values for input from which
  to pick a mean value for each id (applied independently to every Input
  dimension)

- diff_group:

  A number, indicating the mean difference between consecutive groups.

- var_sample:

  A number, indicating the noise variance for each new sample of a id

## Value

A full dataset of synthetic data, with columns `ID`, `Group`, `Sample`,
`Input_ID`, `Input`, `Output`.

## Examples

``` r
data <- simu_db()
```
