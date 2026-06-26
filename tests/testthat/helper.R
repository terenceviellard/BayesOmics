library(methods)

# Prevent plotting tests from opening a real graphics device (Rplots.pdf)
# during devtools::test()/R CMD check.
grDevices::pdf(NULL)

# ── Kernel ──────────────────────────────────────────────────────────────────

make_kernel <- function(hp = c(1.0, 1.0)) {
  k <- new("SEKernel")
  keRnel::set_hyperparameters(k, hp)
}

# ── Data ─────────────────────────────────────────────────────────────────────

make_data <- function(nb_id = 5, nb_group = 2, nb_sample = 1, seed = 42) {
  set.seed(seed)
  simu_db(nb_id = nb_id, nb_group = nb_group, nb_sample = nb_sample)
}

# ── Posteriors ───────────────────────────────────────────────────────────────

make_posteriors <- function(nb_id = 5, nb_group = 2, seed = 42) {
  data <- make_data(nb_id = nb_id, nb_group = nb_group, seed = seed)
  multi_posterior_mean(data, make_kernel())
}

# ── sample_distrib (for plot_distrib) ────────────────────────────────────────

make_sample_distrib <- function(groups = c("G1", "G2"),
                                ids    = c("ID_1", "ID_2"),
                                n      = 200,
                                seed   = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(groups), function(gi) {
    g <- groups[gi]
    do.call(rbind, lapply(ids, function(id) {
      data.frame(
        ID      = id,
        Group   = g,
        Sample  = stats::rnorm(n, mean = gi * 3),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# ── Hand-crafted posterior list (multi_posterior_mean() output format) ───────
# Format: list(kernels = <named list of matrices keyed by Input value strings>,
#              groups  = <named list per group with muk, id_to_input, kernel_key, scale>)

make_simple_results <- function(n_groups = 2, n_ids = 2) {
  groups   <- paste0("g", seq_len(n_groups))
  ids      <- paste0("ID_", seq_len(n_ids))
  inputs   <- seq_len(n_ids) * 1.0
  kern_mat <- diag(n_ids)
  dimnames(kern_mat) <- list(
    BayesOmics:::format_input_key(inputs),
    BayesOmics:::format_input_key(inputs)
  )

  list(
    kernels = list(k = kern_mat),
    groups  = stats::setNames(
      lapply(seq_len(n_groups), function(g) {
        list(
          muk         = stats::setNames(seq_len(n_ids) * g * 1.0, ids),
          id_to_input = stats::setNames(inputs, ids),
          kernel_key  = "k",
          scale       = 1
        )
      }),
      groups
    )
  )
}

# ── Custom posterior list from explicit (muk, sigma) pairs per group ─────────
# groups_spec: named list, each element = list(
#   muk   = named numeric vector,
#   sigma = the RAW kernel matrix (rows/cols in the same order as names(muk)),
#   scale = optional divisor applied to `sigma` to get the posterior covariance
#           (defaults to 1, i.e. sigma IS the posterior covariance).
# )
# Groups whose `sigma` matrices are numerically equal are assigned the SAME
# kernel_key (mirroring multi_posterior_mean()'s deduplication-by-Input-set),
# so that scale-only differences between groups exercise the c != 1 branch of
# calculate_group_overlaps() while still sharing one raw kernel matrix. Groups
# with genuinely different `sigma` values get distinct keys, which is what
# triggers calculate_group_overlaps()'s "different kernel_key" error.

make_custom_results <- function(groups_spec) {
  kernels <- list()
  groups  <- list()
  for (gname in names(groups_spec)) {
    spec   <- groups_spec[[gname]]
    ids    <- names(spec$muk)
    inputs <- seq_along(ids) * 1.0
    raw    <- spec$sigma
    scale  <- if (!is.null(spec$scale)) spec$scale else 1

    key <- NULL
    for (existing_key in names(kernels)) {
      existing <- kernels[[existing_key]]
      if (identical(dim(existing), dim(raw)) &&
          isTRUE(all.equal(unname(existing), unname(raw)))) {
        key <- existing_key
        break
      }
    }
    if (is.null(key)) {
      key  <- paste0("k", length(kernels) + 1)
      kmat <- raw
      dimnames(kmat) <- list(
        BayesOmics:::format_input_key(inputs),
        BayesOmics:::format_input_key(inputs)
      )
      kernels[[key]] <- kmat
    }

    groups[[gname]] <- list(
      muk         = spec$muk,
      id_to_input = stats::setNames(inputs, ids),
      kernel_key  = key,
      scale       = scale
    )
  }
  list(kernels = kernels, groups = groups)
}
