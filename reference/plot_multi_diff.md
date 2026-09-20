# Plot a multivariate summary of group differences

For every pair of groups present in `multi_diff`, plots the empirical
distribution of the number of IDs for which group1's posterior draw
exceeds group2's (see
[`compute_multi_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_multi_diff.md)),
arranged as an upper-triangular grid of panels via
[`gridExtra::grid.arrange`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
– one panel per pair, with no cap on the number of groups (mirrors
[`plot_distrib_each_pair`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib_each_pair.md)'s
equivalent no-cap behaviour). This complements
[`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md)'s
per-id view with a region-wide, uncertainty-aware summary of whether two
groups are differential.

## Usage

``` r
plot_multi_diff(multi_diff, plot_mean = TRUE, cumulative = FALSE)
```

## Arguments

- multi_diff:

  A list, typically coming from
  [`compute_multi_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_multi_diff.md),
  with elements `Diff_proba` and `Diff_mean` (and optionally
  `Overlap_coef`).

- plot_mean:

  A boolean, indicating whether an additional panel showing the
  posterior mean of every ID (coloured by group, see
  [`plot_posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/plot_posterior_mean.md))
  should be displayed. Defaults to `TRUE`.

- cumulative:

  A boolean, indicating whether each panel shows the cumulative
  distribution (`TRUE`) or the probability mass (`FALSE`, default) of
  the number of IDs where group1's posterior draw exceeds group2's.

## Value

The result of
[`gridExtra::grid.arrange`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html):
a grid of panels, one per group pair, plus (if `plot_mean = TRUE`) the
id-mean panel.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 3, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
multi_diff <- compute_multi_diff(samples, results = posterior)
plot_multi_diff(multi_diff)
```
