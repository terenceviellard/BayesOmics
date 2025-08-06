#' @importFrom methods is new show slot<- slotNames
#' @importFrom stats runif
NULL


# Utils -------------------------------------------------------------------
# Fonction pour s'assurer que l'objet est un AbstractKernel
ensure_abstract_kernel <- function(kernel) {
  if (!is(kernel, "AbstractKernel")) {
    return(new("ConstantKernel", value = kernel))
  }
  return(kernel)
}

# Fonction pour obtenir les noms des hyperparamètres
get_hyperparameter_names <- function(kernel) {
  hps <- gt_HPs(kernel)
  names(unlist(hps))
}

# Fonction pour obtenir les valeurs des hyperparamètres
get_hyperparameter_values <- function(kernel) {
  hps <- gt_HPs(kernel)
  unlist(hps)
}


#' @title Replace Hps Kernels
#'
#' @description
#' Update HPs of complex kernels
#'
#' @param kernel Previous kernel
#' @param values New HPs
#'
#' @return A kernel with new HPs
#' @export
#'
set_hyperparameters <- function(kernel, values) {
  slot_names <- get_hyperparameter_names(kernel)

  # Check if the values are named
  if (is.null(names(values))) {
    if (length(values) != length(slot_names)) {
      stop("The number of values provided does not match the number of slots in the kernel.")
    }
    names(values) <- slot_names
  } else {
    for (name in names(values)) {
      if (!(name %in% slot_names)) {
        stop(paste("The slot", name, "does not exist in the kernel."))
      }
    }
  }

  # Update the hyperparameters
  for (name in names(values)) {
    value <- values[[name]]
    if (is(kernel, "SumKernel") || is(kernel, "ProductKernel")) {
      for (i in seq_along(kernel@kernels)) {
        k <- kernel@kernels[[i]]
        if (name %in% get_hyperparameter_names(k)) {
          # Directly set the slot value if the hyperparameter name exists
          slot(k, name) <- value

          # Update the kernel in the composite kernel
          kernel@kernels[[i]] <- k

        }
      }
    } else {
      if (name %in% slotNames(kernel)) {
        slot(kernel, name) <- value
      }
    }
  }


  return(kernel)
}


# AbstractKernel ----------------------------------------------------------
setClass("AbstractKernel",
         representation = representation(),
         prototype = prototype())

setGeneric("pairwise_kernel", function(obj, x, y) standardGeneric("pairwise_kernel"))
setGeneric("pretty_print", function(obj) standardGeneric("pretty_print"))
setGeneric("gt_HPs", function(obj) standardGeneric("gt_HPs"))
setGeneric("kernel_deriv", function(obj, x, y, param) standardGeneric("kernel_deriv"))
setGeneric("kernel_deriv_exp", function(obj, x, y, param) standardGeneric("kernel_deriv_exp"))

setMethod("pairwise_kernel", "AbstractKernel",
          function(obj, x, y) {
            stop("pairwise_kernel method must be implemented in subclass")
          })

setMethod("show", "AbstractKernel", function(object) {
  cat("Abstract Kernel Object\n")
})

setMethod("pretty_print", "AbstractKernel", function(obj) {
  "AbstractKernel()"
})

setMethod("gt_HPs", "AbstractKernel",
          function(obj) {
            stop("gt_HPs method must be implemented in subclass")
          })

setMethod("kernel_deriv", "AbstractKernel",
          function(obj, x, y, param) {
            stop("kernel_deriv method must be implemented in subclass")
          })

# SumKernel ---------------------------------------------------------------
setClass("SumKernel",
         contains = "AbstractKernel",
         slots = c(kernels = "list"))

setMethod("initialize", "SumKernel",
          function(.Object, kernels) {
            .Object@kernels <- lapply(kernels, ensure_abstract_kernel)
            return(.Object)
          })

setMethod("pairwise_kernel", "SumKernel",
          function(obj, x, y) {
            result <- 0
            for (kernel in obj@kernels) {
              result <- result + pairwise_kernel(kernel, x, y)
            }
            return(result)
          })

setMethod("show", "SumKernel", function(object) {
  cat("Sum Kernel:\n")
  for (i in seq_along(object@kernels)) {
    cat("  Kernel", i, ":\n")
    show(object@kernels[[i]])
  }
})

setMethod("pretty_print", "SumKernel", function(obj) {
  kernel_strings <- sapply(obj@kernels, pretty_print)
  paste0("[", paste(kernel_strings, collapse = " + "), "]")
})

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 ANY.
#' @return An object of class AbstractKernel.
#' @export
setMethod("+", signature(e1 = "AbstractKernel", e2 = "ANY"),
          function(e1, e2) {
            new("SumKernel", kernels = list(e1, ensure_abstract_kernel(e2)))
          })

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 ANY.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class AbstractKernel.
#' @export
setMethod("+", signature(e1 = "ANY", e2 = "AbstractKernel"),
          function(e1, e2) {
            new("SumKernel", kernels = list(ensure_abstract_kernel(e1), e2))
          })

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class AbstractKernel.
#' @export
setMethod("+", signature(e1 = "AbstractKernel", e2 = "AbstractKernel"),
          function(e1, e2) {
            new("SumKernel", kernels = list(e1, e2))
          })

setMethod("gt_HPs", "SumKernel",
          function(obj) {
            hps <- list()
            for (kernel in obj@kernels) {
              hps <- c(hps, gt_HPs(kernel))
            }
            return(hps)
          })

setMethod("kernel_deriv", "SumKernel",
          function(obj, x, y, param) {
            deriv <- matrix(0, nrow = nrow(x), ncol = nrow(y))
            for (kernel in obj@kernels) {
              if (param %in% names(gt_HPs(kernel))) {
                deriv <- deriv + kernel_deriv(kernel, x, y, param)
              }
            }
            return(deriv)
          })

# ProductKernel -----------------------------------------------------------
setClass("ProductKernel",
         contains = "AbstractKernel",
         slots = c(kernels = "list"))

setMethod("initialize", "ProductKernel",
          function(.Object, kernels) {
            .Object@kernels <- lapply(kernels, ensure_abstract_kernel)
            return(.Object)
          })

setMethod("pairwise_kernel", "ProductKernel",
          function(obj, x, y) {
            result <- 1
            for (kernel in obj@kernels) {
              result <- result * pairwise_kernel(kernel, x, y)
            }
            return(result)
          })

setMethod("show", "ProductKernel", function(object) {
  cat("Product Kernel:\n")
  for (i in seq_along(object@kernels)) {
    cat("  Kernel", i, ":\n")
    show(object@kernels[[i]])
  }
})

setMethod("pretty_print", "ProductKernel", function(obj) {
  kernel_strings <- sapply(obj@kernels, pretty_print)
  paste0("[", paste(kernel_strings, collapse = " * "), "]")
})

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 ANY.
#' @return An object of class AbstractKernel.
#' @export
setMethod("*", signature(e1 = "AbstractKernel", e2 = "ANY"),
          function(e1, e2) {
            new("ProductKernel", kernels = list(e1, ensure_abstract_kernel(e2)))
          })

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 ANY.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class AbstractKernel.
#' @export
setMethod("*", signature(e1 = "ANY", e2 = "AbstractKernel"),
          function(e1, e2) {
            new("ProductKernel", kernels = list(ensure_abstract_kernel(e1), e2))
          })

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class AbstractKernel.
#' @export
setMethod("*", signature(e1 = "AbstractKernel", e2 = "AbstractKernel"),
          function(e1, e2) {
            new("ProductKernel", kernels = list(e1, e2))
          })

setMethod("gt_HPs", "ProductKernel",
          function(obj) {
            hps <- list()
            for (kernel in obj@kernels) {
              hps <- c(hps, gt_HPs(kernel))
            }
            return(hps)
          })

setMethod("kernel_deriv", "ProductKernel",
          function(obj, x, y, param) {
            deriv <- matrix(0, nrow = nrow(x), ncol = nrow(y))
            for (i in seq_along(obj@kernels)) {
              kernel <- obj@kernels[[i]]
              if (param %in% names(gt_HPs(kernel))) {
                kernel_deriv_value <- kernel_deriv(kernel, x, y, param)
                other_kernels_value <- 1
                for (j in seq_along(obj@kernels)) {
                  if (i != j) {
                    other_kernels_value <- other_kernels_value * pairwise_kernel(obj@kernels[[j]], x, y)
                  }
                }
                deriv <- deriv + other_kernels_value * kernel_deriv_value
              }
            }
            return(deriv)
          })

# ConstantKernel ----------------------------------------------------------
setClass("ConstantKernel",
         contains = "AbstractKernel",
         slots = c(value_c = "numeric"))

setMethod("initialize", "ConstantKernel",
          function(.Object, value_c = runif(1, 0, 3)) {
            .Object@value_c <- value_c
            return(.Object)
          })

setMethod("pairwise_kernel", "ConstantKernel",
          function(obj, x, y) {
            matrix(obj@value_c, nrow = nrow(x), ncol = nrow(y))
          })

setMethod("show", "ConstantKernel", function(object) {
  cat("Constant Kernel:\n")
  cat("  Value:", object@value_c, "\n")
})

setMethod("pretty_print", "ConstantKernel", function(obj) {
  sprintf("ConstantKernel(%.2f)", obj@value_c)
})

setMethod("gt_HPs", "ConstantKernel",
          function(obj) {
            list(value_c = obj@value_c)
          })

setMethod("kernel_deriv", "ConstantKernel",
          function(obj, x, y, param) {
            if (param == "value_c") {
              return(matrix(1, nrow = nrow(x), ncol = nrow(y)))
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

# NoiseKernel ------------------------------------------------------------
setClass("NoiseKernel",
         contains = "AbstractKernel",
         slots = c(value_c = "numeric"))

setMethod("initialize", "NoiseKernel",
          function(.Object, value_c = runif(1, 0, 3)) {
            .Object@value_c <- value_c
            return(.Object)
          })

setMethod("pairwise_kernel", "NoiseKernel",
          function(obj, x, y) {
            obj@value_c * diag(nrow(x))
          })

setMethod("show", "NoiseKernel", function(object) {
  cat("Noise Kernel:\n")
  cat("  Value:", object@value_c, "\n")
})

setMethod("pretty_print", "NoiseKernel", function(obj) {
  sprintf("NoiseKernel(%.2f)", obj@value_c)
})

setMethod("gt_HPs", "NoiseKernel",
          function(obj) {
            list(value_c = obj@value_c)
          })

setMethod("kernel_deriv", "NoiseKernel",
          function(obj, x, y, param) {
            if (param == "value_c") {
              return(diag(nrow(x)))
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

# SEKernel ----------------------------------------------------------------
setClass("SEKernel",
         contains = "AbstractKernel",
         slots = c(variance_se = "numeric", length_scale_se = "numeric"))

setMethod("initialize", "SEKernel",
          function(.Object, variance_se = runif(1, 0, 3), length_scale_se = runif(1, 0, 3)) {
            .Object@variance_se <- variance_se
            .Object@length_scale_se <- length_scale_se
            return(.Object)
          })

setMethod("pairwise_kernel", "SEKernel",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            return(obj@variance_se * exp(-dx / (2 * obj@length_scale_se^2)))
          })

setMethod("kernel_deriv", "SEKernel",
          function(obj, x, y, param) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            if (param == "variance_se") {
              return(pairwise_kernel(obj, x, y) / obj@variance_se)
            } else if (param == "length_scale_se") {
              return(obj@variance_se * exp(-dx / (2 * obj@length_scale_se^2)) * dx / (obj@length_scale_se^3))
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "SEKernel", function(object) {
  cat("Squared Exponential Kernel:\n")
  cat("  Variance:", object@variance_se, "\n")
  cat("  Length Scale:", object@length_scale_se, "\n")
})

setMethod("pretty_print", "SEKernel", function(obj) {
  sprintf("SEKernel(variance=%.2f, length_scale=%.2f)", obj@variance_se, obj@length_scale_se)
})

setMethod("gt_HPs", "SEKernel",
          function(obj) {
            list(variance_se = obj@variance_se, length_scale_se = obj@length_scale_se)
          })

# LinearKernel -------------------------------------------------------------
setClass("LinearKernel",
         contains = "AbstractKernel",
         slots = c(sigma2_b = "numeric", sigma2_v = "numeric", c = "numeric"))

setMethod("initialize", "LinearKernel",
          function(.Object, sigma2_b = runif(1, 0, 3), sigma2_v = runif(1, 0, 3), c = runif(1, 0, 3)) {
            .Object@sigma2_b <- sigma2_b
            .Object@sigma2_v <- sigma2_v
            .Object@c <- c
            return(.Object)
          })

setMethod("pairwise_kernel", "LinearKernel",
          function(obj, x, y) {
            x_centered <- x - obj@c
            y_centered <- y - obj@c
            product <- tcrossprod(x_centered, y_centered)
            return(obj@sigma2_b + obj@sigma2_v * product)
          })

setMethod("kernel_deriv", "LinearKernel",
          function(obj, x, y, param) {
            if (param == "sigma2_b") {
              return(matrix(1, nrow = nrow(x), ncol = nrow(y)))
            } else if (param == "sigma2_v") {
              x_centered <- x - obj@c
              y_centered <- y - obj@c
              return(tcrossprod(x_centered, y_centered))
            } else if (param == "c") {
              x <- -obj@sigma2_v * (outer(x, y, FUN = "+") - 2 * obj@c)
              return(x[, 1, , drop = TRUE])
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "LinearKernel", function(object) {
  cat("Linear Kernel:\n")
  cat("  Sigma squared b:", object@sigma2_b, "\n")
  cat("  Sigma squared v:", object@sigma2_v, "\n")
  cat("  c:", object@c, "\n")
})

setMethod("pretty_print", "LinearKernel", function(obj) {
  sprintf("LinearKernel(sigma2_b=%.2f, sigma2_v=%.2f, c=%.2f)", obj@sigma2_b, obj@sigma2_v, obj@c)
})

setMethod("gt_HPs", "LinearKernel",
          function(obj) {
            list(sigma2_b = obj@sigma2_b, sigma2_v = obj@sigma2_v, c = obj@c)
          })

# RationalQuadraticKernel -------------------------------------------------
setClass("RationalQuadraticKernel",
         contains = "AbstractKernel",
         slots = c(variance_rq = "numeric", length_scale_rq = "numeric", alpha_rq = "numeric"))

setMethod("initialize", "RationalQuadraticKernel",
          function(.Object, variance_rq = runif(1, 0, 3), length_scale_rq = runif(1, 0, 3), alpha_rq = runif(1, 0, 3)) {
            .Object@variance_rq <- variance_rq
            .Object@length_scale_rq <- length_scale_rq
            .Object@alpha_rq <- alpha_rq
            return(.Object)
          })

setMethod("pairwise_kernel", "RationalQuadraticKernel",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            return(obj@variance_rq * (1 + dx / (2 * obj@alpha_rq * obj@length_scale_rq^2))^(-obj@alpha_rq))
          })

setMethod("kernel_deriv", "RationalQuadraticKernel",
          function(obj, x, y, param) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            if (param == "variance_rq") {
              return(pairwise_kernel(obj, x, y) / obj@variance_rq)
            } else if (param == "length_scale_rq") {
              term <- pairwise_kernel(obj, x, y) * dx / ((obj@length_scale_rq^3) * (1 + dx / (2 * obj@alpha_rq * obj@length_scale_rq^2)))
              return(term)
            } else if (param == "alpha_rq") {
              term1 <- (2 * obj@alpha_rq * obj@length_scale_rq^2 + dx)
              term2 <- log(1 + dx / (2 * obj@alpha_rq * obj@length_scale_rq^2))
              term3 <- dx / term1
              combined_term <- term2 - term3
              kernel_value <- pairwise_kernel(obj, x, y)
              deriv <- -kernel_value * combined_term
              return(deriv)
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "RationalQuadraticKernel", function(object) {
  cat("Rational Quadratic Kernel:\n")
  cat("  Variance:", object@variance_rq, "\n")
  cat("  Length Scale:", object@length_scale_rq, "\n")
  cat("  Alpha:", object@alpha_rq, "\n")
})

setMethod("pretty_print", "RationalQuadraticKernel", function(obj) {
  sprintf("RationalQuadraticKernel(variance=%.2f, length_scale=%.2f, alpha=%.2f)", obj@variance_rq, obj@length_scale_rq, obj@alpha_rq)
})

setMethod("gt_HPs", "RationalQuadraticKernel",
          function(obj) {
            list(variance_rq = obj@variance_rq, length_scale_rq = obj@length_scale_rq, alpha_rq = obj@alpha_rq)
          })

# PeriodicKernel ----------------------------------------------------------
setClass("PeriodicKernel",
         contains = "AbstractKernel",
         slots = c(variance_per = "numeric", length_scale_per = "numeric", period = "numeric"))

setMethod("initialize", "PeriodicKernel",
          function(.Object, variance_per = runif(1, 0, 3), length_scale_per = runif(1, 0, 3), period = runif(1, 0, 2 * pi)) {
            .Object@variance_per <- variance_per
            .Object@length_scale_per <- length_scale_per
            .Object@period <- period
            return(.Object)
          })

setMethod("pairwise_kernel", "PeriodicKernel",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            return(obj@variance_per * exp(-2 * (sin(pi * dx / obj@period)^2) / obj@length_scale_per^2))
          })

setMethod("kernel_deriv", "PeriodicKernel",
          function(obj, x, y, param) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            if (param == "variance_per") {
              return(exp(-2 * (sin(pi * dx / obj@period)^2) / obj@length_scale_per^2) / obj@variance_per)
            } else if (param == "length_scale_per") {
              x <- 4 * obj@variance_per * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 * exp(-2 * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 / obj@length_scale_per^2) / obj@length_scale_per^3
              return(x[, 1, , drop = TRUE])
            } else if (param == "period") {
              x <- 2 * pi * obj@variance_per * outer(x, y, FUN = "-") * sin(2 * pi * outer(x, y, FUN = "-") / obj@period) * exp(-2 * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 / obj@length_scale_per^2) / (obj@length_scale_per^2 * obj@period^2)
              return(x[, 1, , drop = TRUE])
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "PeriodicKernel", function(object) {
  cat("Periodic Kernel:\n")
  cat("  Variance:", object@variance_per, "\n")
  cat("  Length Scale:", object@length_scale_per, "\n")
  cat("  Period:", object@period, "\n")
})

setMethod("pretty_print", "PeriodicKernel", function(obj) {
  sprintf("PeriodicKernel(variance=%.2f, length_scale=%.2f, period=%.2f)", obj@variance_per, obj@length_scale_per, obj@period)
})

setMethod("gt_HPs", "PeriodicKernel",
          function(obj) {
            list(variance_per = obj@variance_per, length_scale_per = obj@length_scale_per, period = obj@period)
          })

# MaternKernel 1/2 ------------------------------------------------------------
setClass("MaternKernel12",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric"))

setMethod("initialize", "MaternKernel12",
          function(.Object, length_scale_mat = runif(1, 0, 3)) {
            .Object@length_scale_mat <- length_scale_mat
            return(.Object)
          })

setMethod("pairwise_kernel", "MaternKernel12",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            return(exp(-dx / obj@length_scale_mat))
          })

setMethod("kernel_deriv", "MaternKernel12",
          function(obj, x, y, param) {
            if (param == "length_scale_mat") {
              dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
              return((dx / obj@length_scale_mat^2) * exp(-dx / obj@length_scale_mat))
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "MaternKernel12", function(object) {
  cat("MaternKernel 1/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

setMethod("pretty_print", "MaternKernel12", function(obj) {
  sprintf("MaternKernel12(length_scale=%.2f)", obj@length_scale_mat)
})

setMethod("gt_HPs", "MaternKernel12",
          function(obj) {
            list(length_scale_mat = obj@length_scale_mat)
          })

# MaternKernel 3/2 ------------------------------------------------------------
setClass("MaternKernel32",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric"))

setMethod("initialize", "MaternKernel32",
          function(.Object, length_scale_mat = runif(1, 0, 3)) {
            .Object@length_scale_mat <- length_scale_mat
            return(.Object)
          })

setMethod("pairwise_kernel", "MaternKernel32",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            sqrt3_r_div_l <- (sqrt(3) * dx) / obj@length_scale_mat
            return((1 + sqrt3_r_div_l) * exp(-sqrt3_r_div_l))
          })

setMethod("kernel_deriv", "MaternKernel32",
          function(obj, x, y, param) {
            if (param == "length_scale_mat") {
              dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
              return(3 * dx^2 * exp(-sqrt(3) * dx / obj@length_scale_mat) / obj@length_scale_mat^3)
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "MaternKernel32", function(object) {
  cat("MaternKernel 3/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

setMethod("pretty_print", "MaternKernel32", function(obj) {
  sprintf("MaternKernel32(length_scale=%.2f)", obj@length_scale_mat)
})

setMethod("gt_HPs", "MaternKernel32",
          function(obj) {
            list(length_scale_mat = obj@length_scale_mat)
          })

# MaternKernel 5/2 ------------------------------------------------------------
setClass("MaternKernel52",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric"))

setMethod("initialize", "MaternKernel52",
          function(.Object, length_scale_mat = runif(1, 0, 3)) {
            .Object@length_scale_mat <- length_scale_mat
            return(.Object)
          })

setMethod("pairwise_kernel", "MaternKernel52",
          function(obj, x, y) {
            dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
            sqrt5_r_div_l <- (sqrt(5) * dx) / obj@length_scale_mat
            return((1 + sqrt5_r_div_l + (5.0 / 3.0) * (dx / obj@length_scale_mat)^2) * exp(-sqrt5_r_div_l))
          })

setMethod("kernel_deriv", "MaternKernel52",
          function(obj, x, y, param) {
            if (param == "length_scale_mat") {
              dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
              return(5 * dx^2 * exp(-sqrt(5) * dx / obj@length_scale_mat) * (obj@length_scale_mat + sqrt(5) * dx)^3 / obj@length_scale_mat^4)
            } else {
              stop("Unknown parameter for derivative calculation.")
            }
          })

setMethod("show", "MaternKernel52", function(object) {
  cat("MaternKernel 5/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

setMethod("pretty_print", "MaternKernel52", function(obj) {
  sprintf("MaternKernel52(length_scale=%.2f)", obj@length_scale_mat)
})

setMethod("gt_HPs", "MaternKernel52",
          function(obj) {
            list(length_scale_mat = obj@length_scale_mat)
          })




# OPITM_HP ----------------------------------------------------------------



#' @noRd
chol_inv_jitter <- function(mat, pen_diag) {

  diag(mat) <- diag(mat) + pen_diag

  tryCatch(
    {
      chol_mat <- chol(mat)
      inv_mat <- chol2inv(chol_mat)
      inv_mat
    },
    error = function(e) {
      chol_inv_jitter(mat, 10 * pen_diag)
    }
  )
}


#' @noRd
dmnorm <- function(x, mu, inv_Sigma, log = FALSE) {
  # Ensure x is a matrix
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  n <- nrow(x)
  p <- ncol(x)

  # Ensure mu is compatible with x
  if (is.vector(mu)) {
    if (length(mu) != p) {
      stop("The length of mu must match the number of columns in x.")
    }
    mu <- matrix(rep(mu, n), ncol = p, byrow = TRUE)
  } else if (is.matrix(mu)) {
    if (ncol(mu) != p || nrow(mu) != n) {
      stop("The dimensions of mu must match the dimensions of x.")
    }
  } else {
    stop("mu must be a vector or a matrix.")
  }

  # Compute z as the difference between x and mu
  z <- x - mu


  # Ensure dimensions are compatible for multiplication
  if (ncol(z) != nrow(inv_Sigma)) {
    if (nrow(z) == nrow(inv_Sigma)) {
      # Transpose z if necessary
      z <- t(z)
    } else {
      stop("The number of columns in z must match the number of rows in inv_Sigma.")
    }
  }

  logdetS <- try(-determinant(inv_Sigma, logarithm = TRUE)$modulus, silent = TRUE)
  attributes(logdetS) <- NULL

  # Perform the multiplication correctly
  ssq <- sum((z %*% inv_Sigma) * z)

  loglik <- as.vector(-(n * log(2 * pi) + logdetS + ssq) / 2)

  if (log) {
    return(loglik)
  } else {
    return(exp(loglik))
  }
}






#' @noRd
sum_logGaussian <- function(hp, db, mean, kern, post_cov, pen_diag) {
  kern <- set_hyperparameters(kern, hp)
  input <- db$Input
  input <- as.matrix(input)
  cov <- pairwise_kernel(kern, input, input) + post_cov
  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  # Ensure db$Output is a matrix or data frame
  if (!is.matrix(db$Output) && !is.data.frame(db$Output)) {
    # If db$Output is a vector, convert it to a matrix with a single column
    if (is.vector(db$Output)) {
      db$Output <- matrix(db$Output, ncol = 1)
    } else {
      stop("db$Output must be a matrix, data frame, or vector.")
    }
  }

  ncol_output <- ncol(db$Output)
  if (ncol_output == 0) {
    stop("db$Output must have at least one column.")
  }

  # Ensure mean is compatible with the number of columns in db$Output
  if (is.vector(mean)) {
    if (length(mean) == 1) {
      mean <- rep(mean, times = ncol_output)
    } else if (length(mean) != ncol_output) {
      stop("The length of mean must match the number of columns in db$Output.")
    }
  } else if (is.matrix(mean)) {
    if (ncol(mean) != ncol_output || nrow(mean) != nrow(db$Output)) {
      stop("The dimensions of mean must match the dimensions of db$Output.")
    }
  } else {
    stop("mean must be a vector or a matrix.")
  }

  log_likelihoods <- dmnorm(db$Output, mean, inv, log = TRUE)
  neg_sum_log_likelihood <- -sum(log_likelihoods)
  return(neg_sum_log_likelihood)
}




#' @noRd
gr_sum_logGaussian <- function(hp, db, mean, kern, post_cov, pen_diag) {
  kern <- set_hyperparameters(kern, hp)
  list_hp <- get_hyperparameter_names(kern)

  output <- db$Output
  input <- db$Input
  input <- as.matrix(input)

  cov <- pairwise_kernel(kern, input, input) + post_cov

  inv <- chol_inv_jitter(cov, pen_diag = pen_diag)

  prod_inv <- inv %*% (output - mean)
  common_term <- prod_inv %*% t(prod_inv) - inv

  floop <- function(hp_name) {
    kern_deriv <- kernel_deriv(kern, input, input, hp_name)
    grad_term <- sum(diag((-0.5 * (common_term %*% kern_deriv))))
  }

  sapply(list_hp, floop)
}

#' Optimize Hyperparameters for a kernel
#'
#' This function optimizes hyperparameters using the L-BFGS-B optimization method.
#'
#' @param hp A vector of initial hyperparameters to be optimized.
#' @param db The dataset used for optimization.
#' @param mean The mean function
#' @param kern The kernel function
#' @param post_cov The posterior covariance function.
#' @param pen_diag A penalty term added to the diagonal of the covariance matrix for numerical stability.
#' @param verbose A logical value indicating whether to return the full optimization result or just the optimized parameters. Defaults to FALSE.
#'
#' @return If `verbose` is FALSE, a vector of optimized hyperparameters; otherwise, the full result from the `optim` function.
#' @export
optim_hp <- function(hp, db, mean, kern, post_cov, pen_diag, verbose = FALSE) {
  # Call the optimization function
  result <- stats::optim(
    par = hp,
    fn = function(hp) sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag),
    gr = function(hp) gr_sum_logGaussian(hp, db, mean, kern, post_cov, pen_diag),
    method = "L-BFGS-B",
    control = list(factr = 1e7, maxit = 1000)
  )
  if (!verbose) {
    return(result$par)
  } else {
    return(result)

  }

}
