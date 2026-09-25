# Generate a Synthetic Dataset with Kernel-Structured Covariance

Simulate a complete training dataset, similar to
[`simu_db()`](https://terenceviellard.github.io/BayesOmics/reference/simu_db.md),
but consistent with the generative model underlying
[`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md):
for each group, the vector of Output values across the `nb_id` ids
(indexed by their Input position) is drawn jointly from a multivariate
normal distribution whose covariance is given by the kernel applied
pairwise to the (shared, per-id) Input values, plus independent
measurement noise. Each replicate is an independent draw of that same
distribution: \$\$y_n \mid \mu_g \sim \mathcal{N}(\mu_g,\\
\Sigma\_\theta + \sigma^2 I), \quad n = 1,\dots,N_g\$\$ where
\\\Sigma\_\theta\\ is the kernel matrix over the (shared) Input values
and \\\sigma^2\\ is `var_sample`. This makes \\\Sigma\_\theta\\ – not
just the Input values – shared across every group, exactly as assumed by
[`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md)
and required by
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md).

By default (`mu_random = FALSE`) the per-group mean \\\mu_g\\ is a
deterministic per-id baseline shifted by the corresponding element of
`diff_group`. Setting `mu_random = TRUE` instead draws \\\mu_g\\ from
the conjugate prior, \\\mu_g \sim \mathcal{N}(\mu_0,\\
\Sigma\_\theta/\lambda_0)\\ (plus the group shift) – useful for
validating the posterior's frequentist coverage against the true
\\\mu_g\\.

## Usage

``` r
simu_db_kernel(
  nb_id = 5,
  nb_group = 2,
  nb_sample = 5,
  nb_dim = 1,
  range_output = c(0, 50),
  range_input = c(0, 50),
  diff_group = 3,
  var_sample = 2,
  kernel = NULL,
  mu_random = FALSE,
  mu_0 = NULL,
  lambda_0 = 1,
  pen_diag = 1e-06,
  group_labels = NULL,
  integer_input = FALSE,
  input_grid = FALSE
)
```

## Arguments

- nb_id:

  An integer, indicating the number of features (ids) in the data.

- nb_group:

  An integer, indicating the number of groups/conditions.

- nb_sample:

  An integer or a vector of length `nb_group`. When a scalar, all groups
  receive the same number of replicates. When a vector, element `g` sets
  the number of replicates for group `g` independently, enabling
  unbalanced designs.

- nb_dim:

  An integer, indicating the number of Input dimensions per id. Defaults
  to `1` (a single scalar Input per id, emitted with `Input_ID = 1`).
  For `nb_dim > 1`, the kernel is evaluated directly on the
  `nb_id x nb_dim` matrix of positions (isotropic by default: a single
  shared `length_scale` across every dimension, unless `kernel` itself
  is an `ard_kernel()` or similar). `range_input` is applied
  independently to every dimension (each dimension's column of positions
  is drawn independently, not as globally-distinct D-dimensional
  tuples).

- range_output:

  A 2-element vector; its midpoint is used as the default baseline mean
  (`mu_0`) when `mu_0` is not supplied.

- range_input:

  A 2-element vector, indicating the range from which Input values are
  drawn for each feature and each dimension. Ignored for the integer
  dimension when `integer_input = TRUE`.

- diff_group:

  A numeric scalar or a vector of length `nb_group`. When a scalar,
  group `g` receives a shift of `diff_group * (g - 1)` (linear
  dose-response, group 1 = reference). When a vector, element `g` is the
  absolute offset applied to group `g`, allowing arbitrary non-linear
  group contrasts.

- var_sample:

  A positive number, the variance of the independent measurement noise
  added on top of the kernel-induced covariance.

- kernel:

  A kernel object (from the keRnel package). Defaults to
  `variance_kernel(variance = 1) * se_kernel(length_scale = 1)` when
  `NULL`.

- mu_random:

  If `TRUE`, draw each group's mean vector from the conjugate prior
  \\\mathcal{N}(\mu_0, \Sigma\_\theta/\lambda_0)\\ instead of a
  deterministic constant baseline. Defaults to `FALSE`.

- mu_0:

  Baseline mean shared by every id (deterministic case) or prior mean
  (`mu_random = TRUE`); defaults to the midpoint of `range_output`.

- lambda_0:

  Prior precision scaling (only used when `mu_random = TRUE`); must be a
  single positive number. Defaults to `1`.

- pen_diag:

  Jitter added to the diagonal of the kernel-induced covariance matrix
  for numerical stability. Defaults to `1e-6`.

- group_labels:

  A character or numeric vector of length `nb_group` providing custom
  labels for the `Group` column. Defaults to `NULL`, which uses `"1"`,
  `"2"`, ..., `"nb_group"`.

- integer_input:

  Logical. If `TRUE`, Input values are distinct integers drawn from
  `[ceiling(range_input[1]), floor(range_input[2])]`. If `FALSE`
  (default), Input values are continuous uniform draws. See also
  `input_grid`.

- input_grid:

  Logical. Only used when `integer_input = TRUE`. If `FALSE` (default),
  integers are sampled uniformly without replacement (random positions).
  If `TRUE`, integers are placed on a regular grid of `nb_id` evenly
  spaced positions across the integer range (deterministic layout).

## Value

A data frame of synthetic data with columns `ID`, `Group`, `Sample`,
`Input_ID`, `Input`, `Output`. Two attributes are attached:

- `mu_true`:

  Named list, one entry per group, each a named vector of true
  per-feature means.

- `base_input`:

  An `nb_id x nb_dim` numeric matrix (`rownames = ID`) giving the Input
  position(s) assigned to each feature.

## Examples

``` r
# Default kernel (variance = 1, length_scale = 1):
data <- simu_db_kernel()

# Custom kernel:
ker <- keRnel::variance_kernel(variance = 2) * keRnel::se_kernel(length_scale = 3)
data <- simu_db_kernel(kernel = ker)

# Unbalanced design: 3 replicates for group 1, 15 for group 2:
data <- simu_db_kernel(nb_sample = c(3, 15))

# Arbitrary group contrasts (non-linear dose-response):
data <- simu_db_kernel(nb_group = 4, diff_group = c(0, 5, 5, 10))

# Custom group labels:
data <- simu_db_kernel(nb_group = 2, group_labels = c("Control", "Treatment"))

# Integer inputs on a regular grid:
data <- simu_db_kernel(nb_id = 10, integer_input = TRUE, input_grid = TRUE,
                       range_input = c(0, 100))

# Multi-dimensional Input (D = 2), isotropic kernel:
data <- simu_db_kernel(nb_dim = 2)
```
