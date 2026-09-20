# Group-Difference Metrics for BayesOmics Posteriors

A `distance_metric` is a lightweight S3 object describing *which*
pairwise divergence/distance to compute between two groups' Gaussian
posteriors, without itself holding any data. Pass one to
[`compute_group_diff`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff.md)
(or, for sample-based metrics,
[`compute_group_diff_samples`](https://terenceviellard.github.io/BayesOmics/reference/compute_group_diff_samples.md))
to get a full group x group matrix of values.

Three S3 generics drive the dispatch, each with a default method so a
new metric only needs to implement `evaluate_metric`:

- `evaluate_metric(metric, mu1, mu2, Sigma1, Sigma2, scale1, scale2, same_kernel, ...)`:

  Computes the metric's value for one ordered pair of groups.
  `mu1`/`mu2` are the groups' posterior means (aligned by ID),
  `Sigma1`/`Sigma2` their posterior covariances (from the internal
  `get_sigmak()` helper), `scale1`/`scale2` the `n_obs + lambda_0`
  divisors from
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md),
  and `same_kernel` whether the two groups share the same `kernel_key`.

- `requires_shared_kernel(metric)`:

  `TRUE` if the metric has no valid formula unless the two groups share
  the same `kernel_key` (only
  [`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)).
  Defaults to `FALSE`.

- `is_symmetric_metric(metric)`:

  `TRUE` if swapping the two groups never changes the value
  ([`kl_metric`](https://terenceviellard.github.io/BayesOmics/reference/kl_metric.md)
  and
  [`mahalanobis_metric`](https://terenceviellard.github.io/BayesOmics/reference/mahalanobis_metric.md)
  are not – the latter is evaluated under \\\Sigma_1^{-1}\\
  specifically, so it differs from its swapped counterpart whenever
  \\\Sigma_1 \neq \Sigma_2\\). Defaults to `TRUE`.

## Usage

``` r
evaluate_metric(
  metric,
  mu1,
  mu2,
  Sigma1,
  Sigma2,
  scale1 = NULL,
  scale2 = NULL,
  same_kernel = FALSE,
  ...
)

requires_shared_kernel(metric)

is_symmetric_metric(metric)
```

## Arguments

- metric:

  A `distance_metric` object.

- mu1, mu2:

  Named numeric vectors, the two groups' posterior means (same IDs).

- Sigma1, Sigma2:

  The two groups' posterior covariance matrices (same ID order as
  `mu1`/`mu2`).

- scale1, scale2:

  The two groups' posterior scale (`n_obs + lambda_0`); `NULL` if
  unavailable.

- same_kernel:

  Whether the two groups share the same `kernel_key`.

- ...:

  Passed through by decorator metrics; unused by concrete metrics.

## Value

`evaluate_metric()` returns a single numeric value.
`requires_shared_kernel()`/`is_symmetric_metric()` return a single
logical.
