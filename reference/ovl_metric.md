# Gaussian Overlapping Coefficient (OVL) Metric

The overlapping coefficient between two groups' posteriors, as used by
[`calculate_group_overlaps`](https://terenceviellard.github.io/BayesOmics/reference/calculate_group_overlaps.md).
Same call in every case: when the two groups share the same `kernel_key`
(so \\\Sigma_2 = c\Sigma_1\\ for a scalar \\c =
\mathrm{scale}\_1/\mathrm{scale}\_2\\), the exact closed-form ratio is
used. Otherwise (e.g. `pooled = FALSE`, or two groups fit with genuinely
different kernels), there is no closed form – the log-likelihood-ratio
boundary between the two Gaussians then has a different weight per axis,
not reducible to one (noncentral) chi-squared variable – so a Monte
Carlo/KDE estimate is used automatically instead (`.ovl_mc_fallback()`):
`n_mc` draws simulated directly from each group's \\N(\mu,\Sigma)\\,
then the two-sample identity \\OVL = E\_{f_1}\[\min(1, f_2/f_1)\] =
E\_{f_2}\[\min(1, f_1/f_2)\]\\ with \\f_1\\/\\f_2\\ estimated by a
diagonal Gaussian KDE. Noisier and increasingly unreliable as the number
of shared IDs grows (the curse of dimensionality for KDE) – only reached
when the exact formula does not apply, so it never affects the pooled
case.

## Usage

``` r
ovl_metric(n_mc = 2000)
```

## Arguments

- n_mc:

  Number of Monte Carlo draws used by the KDE fallback (see above);
  irrelevant (unused) whenever the two groups share a `kernel_key`,
  since the exact formula applies then. Defaults to `2000`.

## Value

A `distance_metric` object.
