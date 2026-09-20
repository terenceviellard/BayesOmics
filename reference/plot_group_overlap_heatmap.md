# Plot a heatmap of pairwise group overlap coefficients

For a given id, compute and display an empirical overlapping coefficient
(the area of intersection of the two groups' posterior density
estimates) between every pair of groups present in `sample_distrib`, as
a group x group heatmap. This gives a single, constant-size overview of
how differentiated every pair of groups is for that id, regardless of
how many groups are present – unlike a full grid of one density panel
per pair, which grows as the number of groups squared and becomes
illegible.

## Usage

``` r
plot_group_overlap_heatmap(sample_distrib, id = NULL, digits = 2)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  function, containing the following columns: `ID`, `Group` and
  `Sample`.

- id:

  A character string, the id for which pairwise overlaps are computed.
  If NULL (default), only the first id appearing in `sample_distrib` is
  used.

- digits:

  Number of decimal digits used when displaying the overlap coefficient
  on each tile. Defaults to `2`.

## Value

A `ggplot` object: a tiled heatmap of groups x groups, filled by the
estimated overlapping coefficient (1 = identical distributions, 0 =
fully separated).

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 4, nb_sample = 5, diff_group = 4)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
plot_group_overlap_heatmap(samples, id = unique(samples$ID)[1])
```
