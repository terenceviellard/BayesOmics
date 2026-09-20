# Compute the Multivariate Distribution of Group Differences

For every pair of groups present in `sample_distrib`, computes the
empirical distribution (over posterior draws) of the number of IDs for
which group1's posterior draw exceeds group2's – a multivariate,
uncertainty-aware complement to the single scalar overlapping
coefficient returned by
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md).

This relies on a documented invariant of
[`sample_posterior`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)'s
output: for a given group, its `n` draws are melted from a single
`n x n_ids` matrix (`Sample = as.vector(mat)`), so filtering
`sample_distrib` by `(Group, ID)` recovers draws in a consistent order
across every ID in that group, without needing an explicit draw index
column. This only holds if `sample_distrib`'s row order has not been
shuffled since
[`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
produced it, and every ID within a group contributes the same number of
draws (checked, and an error raised otherwise; the actual per-draw
pairing itself cannot be independently verified from a melted data frame
with no explicit draw index).

## Usage

``` r
compute_multi_diff(sample_distrib, results = NULL)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  function, containing the following columns: `ID`, `Group` and
  `Sample`.

- results:

  An optional list, typically from
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md),
  with elements `kernels` and `groups`. If supplied,
  [`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md)
  is used to attach an exact `Overlap_coef` to the result. If `NULL`
  (default), `Overlap_coef` is omitted.

## Value

A list with elements:

- `Diff_proba`:

  A tibble with columns `Group1`, `Group2`, `Nb_id` (from 0 to the
  number of shared IDs), `Proba` (the probability mass at that count)
  and `Cumul_proba` (its cumulative sum).

- `Diff_mean`:

  A tibble with columns `ID`, `Group` and `Mean`, the posterior mean of
  each (ID, Group) pair.

- `Overlap_coef`:

  A tibble with columns `Group1`, `Group2` and `Overlap_coef`. Only
  present when `results` is supplied.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 3, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
multi_diff <- compute_multi_diff(samples, results = posterior)
multi_diff$Diff_proba
#> # A tibble: 27 × 5
#>    Group1 Group2 Nb_id Proba Cumul_proba
#>    <chr>  <chr>  <int> <dbl>       <dbl>
#>  1 1      2          0     1           1
#>  2 1      2          1     0           1
#>  3 1      2          2     0           1
#>  4 1      2          3     0           1
#>  5 1      2          4     0           1
#>  6 1      2          5     0           1
#>  7 1      2          6     0           1
#>  8 1      2          7     0           1
#>  9 1      2          8     0           1
#> 10 1      3          0     1           1
#> # ℹ 17 more rows
```
