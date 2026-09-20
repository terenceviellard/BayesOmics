# Plot each pairwise group comparison separately

Return one difference-of-means posterior plot per pair of groups present
in `sample_distrib`, as a named list of `ggplot` objects – with no
summary/selection logic, regardless of how many groups (and therefore
pairs) there are. Use this when you want to inspect or save every
pairwise comparison individually (e.g. in a report with one figure per
pair), as opposed to
[`plot_distrib()`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md)'s
automatic top-N summary view for more than two groups.

## Usage

``` r
plot_distrib_each_pair(
  sample_distrib,
  id = NULL,
  prob_CI = 0.95,
  show_prob = TRUE,
  mean_bar = TRUE
)
```

## Arguments

- sample_distrib:

  A data frame, typically coming from the
  [`sample_posterior()`](https://terenceviellard.github.io/BayesOmics/reference/sample_posterior.md)
  function, containing the following columns: `ID`, `Group` and
  `Sample`.

- id:

  A character string, the id to plot. If NULL (default), only the first
  id appearing in `sample_distrib` is used.

- prob_CI:

  A number, between 0 and 1, the level of the Credible Interval. See
  [`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md).

- show_prob:

  A boolean, whether to display the probability labels. See
  [`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md).

- mean_bar:

  A boolean, whether to display the vertical bar at 0. See
  [`plot_distrib`](https://terenceviellard.github.io/BayesOmics/reference/plot_distrib.md).

## Value

A named list of `ggplot` objects, one per group pair, named
`"<group1>_vs_<group2>"`.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 4, nb_sample = 5, diff_group = 4)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
samples <- sample_posterior(posterior, n = 500)
plots <- plot_distrib_each_pair(samples, id = unique(samples$ID)[1])
names(plots)
#> [1] "1_vs_2" "1_vs_3" "1_vs_4" "2_vs_3" "2_vs_4" "3_vs_4"
```
