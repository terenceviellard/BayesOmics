# Compute Overlapping Coefficient between Groups

Computes a symmetric matrix of pairwise overlapping coefficients (OVL)
between every pair of groups.

For two groups sharing the same underlying kernel matrix (same
`kernel_key`, i.e. the same set of Input values), their posterior
covariances are \\\Sigma_1 = \Sigma/(\lambda_0+N_1)\\ and \\\Sigma_2 =
\Sigma/(\lambda_0+N_2)\\ for a shared raw kernel matrix \\\Sigma\\, so
\\\Sigma_2 = c\\\Sigma_1\\ with \\c =
\mathrm{scale}\_1/\mathrm{scale}\_2\\. Writing \\D^2 =
(\mu_2-\mu_1)^\top\Sigma_1^{-1}(\mu_2-\mu_1)\\ and \\d\\ the dimension
(number of shared IDs), the exact Gaussian overlap is used:

- \\c = 1\\: \$\$OVL = 2 \Phi\\\left(-D/2\right)\$\$

- \\c \ne 1\\, with \\\lambda_1 = D^2/(1-c)^2\\, \\\lambda_2 =
  cD^2/(1-c)^2\\, \\t = (D^2 - d(1-c)\ln c)/(1-c)^2\\: \$\$OVL =
  \begin{cases} F\_{\chi^2_d(\lambda_1)}(ct) + 1 -
  F\_{\chi^2_d(\lambda_2)}(t), & 0\<c\<1 \\
  F\_{\chi^2_d(\lambda_2)}(t) + 1 - F\_{\chi^2_d(\lambda_1)}(ct), & c\>1
  \end{cases}\$\$

Two groups with different `kernel_key` (different sets of Input values)
do not share a common raw kernel matrix \\\Sigma\\, so this closed-form
ratio does not apply; an error is raised in that case.

When \\c\\ is not exactly 1 but very close to it (relative difference
below `1e-6`), the \\c \ne 1\\ formula above becomes numerically
unstable (it divides by \\(1-c)^2\\). In that case a warning is issued
and the \\c = 1\\ formula is used instead, with \\\Sigma_1\\ taken from
whichever group has the smaller scale (i.e. the larger, more
conservative posterior covariance).

Cost: each pair of groups requires one matrix inversion of size \\d
\times d\\ (\\d\\ = number of shared IDs), so the total cost is \\O(G^2
d^3)\\ for \\G\\ groups – e.g. 10 groups with 400 shared IDs already
means about 45 inversions of 400x400 matrices. A warning is issued if
this is likely to be slow (see `max_groups_warn`/`max_dim_warn`).

## Usage

``` r
calculate_group_overlaps(results, max_groups_warn = 50, max_dim_warn = 500)
```

## Arguments

- results:

  A list, typically from
  [`posterior_mean`](https://terenceviellard.github.io/BayesOmics/reference/posterior_mean.md),
  with elements `kernels` and `groups` (one entry per group, each with
  `muk`, `id_to_input`, `kernel_key`, `scale`).

- max_groups_warn:

  Emit a warning about the \\O(G^2 d^3)\\ cost above (see Description)
  if the number of groups exceeds this. Defaults to `50`.

- max_dim_warn:

  Emit the same warning if the number of shared IDs \\d\\ exceeds this.
  Defaults to `500`.

## Value

A symmetric matrix of OVL coefficients in \\\[0, 1\]\\, with 1 on the
diagonal.

## Examples

``` r
data <- simu_db(nb_id = 8, nb_group = 2, nb_sample = 5, diff_group = 5)
kern <- keRnel::variance_kernel(variance = 1) * keRnel::se_kernel(length_scale = 1)
posterior <- posterior_mean(data, kern)
calculate_group_overlaps(posterior)
#>              1            2
#> 1 1.000000e+00 1.267146e-47
#> 2 1.267146e-47 1.000000e+00
```
