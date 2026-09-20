# Total Variation Distance

Total variation distance between the two groups' posteriors, \\TV(f,g) =
1 - \int \min(f,g)\\dx\\. This is a general identity for any two
densities (not an approximation specific to Gaussians): since
\\\int(f+g)\\dx = 2\\, \\\int\|f-g\|\\dx = 2 - 2\int\min(f,g)\\dx\\, so
\\TV = \tfrac12\int\|f-g\|\\dx = 1 - \int\min(f,g)\\dx\\ always. Because
[`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
already computes \\\int\min(f,g)\\dx\\ in closed form, TVD is just
`1 - ` that value – a thin decorator around it rather than a new
closed-form derivation. Forwards
[`requires_shared_kernel`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)/[`is_symmetric_metric`](https://terenceviellard.github.io/BayesOmics/reference/evaluate_metric.md)
to the wrapped metric.

## Usage

``` r
tvd_metric(base = ovl_metric())
```

## Arguments

- base:

  A `distance_metric` object whose value is an overlap coefficient
  ([`ovl_metric`](https://terenceviellard.github.io/BayesOmics/reference/ovl_metric.md)
  by default). Passing anything else is meaningless (TVD is only defined
  via an overlap coefficient), but left as a parameter so a decorated
  OVL variant (e.g. `marginal_metric(ovl_metric())`) can be turned into
  its own TVD analogue too.

## Value

A `distance_metric` object.
