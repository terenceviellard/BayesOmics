#' @importFrom methods is new show slot<- slotNames
#' @importFrom stats runif


# Utils Kernel -------------------------------------------------------------------
#' @title Ensure an Object is an AbstractKernel
#'
#' @description
#' This function checks if the provided object is an instance of the `AbstractKernel` class.
#' If not, it wraps the object in a `ConstantKernel`, provided the object is numeric.
#'
#' @param kernel An object to be checked or wrapped. If not an `AbstractKernel`, it should be a numeric value suitable for initializing a `ConstantKernel`.
#' @return An object of class `AbstractKernel`.
#' @export
#' @examples
#' # Assuming ConstantKernel is a subclass of AbstractKernel
#' kernel <- ensure_abstract_kernel(2.5)
ensure_abstract_kernel <- function(kernel) {
  if (!is(kernel, "AbstractKernel")) {
    if (!is.numeric(kernel)) {
      stop("The provided kernel is not an AbstractKernel and cannot be converted because it is not numeric.")
    }
    return(new("ConstantKernel", value = kernel))
  }
  return(kernel)
}

#' @title Get Hyperparameter Names
#'
#' @description
#' Retrieves the names of hyperparameters from a kernel object.
#'
#' @param kernel An object of class `AbstractKernel` or one of its subclasses.
#' @return A character vector of hyperparameter names.
#' @export
#' @examples
#' # Assuming SEKernel is a subclass of AbstractKernel with hyperparameters
#' kernel <- new("SEKernel", variance_se = 1, length_scale_se = 1)
#' get_hyperparameter_names(kernel)
get_hyperparameter_names <- function(kernel) {
  hps <- gt_HPs(kernel)
  names(unlist(hps))
}

#' @title Get Hyperparameter Values
#'
#' @description
#' Retrieves the values of hyperparameters from a kernel object.
#'
#' @param kernel An object of class `AbstractKernel` or one of its subclasses.
#' @return A named list of hyperparameter values.
#' @export
#' @examples
#' # Assuming SEKernel is a subclass of AbstractKernel with hyperparameters
#' kernel <- new("SEKernel", variance_se = 1, length_scale_se = 1)
#' get_hyperparameter_values(kernel)
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
#' @examples
#' # Assuming SEKernel is a subclass of AbstractKernel with hyperparameters
#' kernel <- new("SEKernel", variance_se = 1, length_scale_se = 1)
#' updated_kernel <- set_hyperparameters(kernel, c(variance_se = 2, length_scale_se = 0.5))
#'
#' @export
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

#' @title Abstract Kernel Class
#' @description
#' This is an abstract class representing a generic kernel. It serves as a base class
#' for specific kernel implementations, providing a common interface for kernel operations.
#'
#'
#' @export
setClass("AbstractKernel",
         representation = representation(),
         prototype = prototype()
)

#' @title Generic Function for Pairwise Kernel Computation
#' @description
#' This generic function computes the pairwise kernel matrix between two sets of vectors.
#'
#' @param obj An object inheriting from the `AbstractKernel` class.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setGeneric("pairwise_kernel", function(obj, x, y) standardGeneric("pairwise_kernel"))

#' @title Generic Function for Pretty Printing
#' @description
#' This generic function provides a pretty-printed string representation of a kernel object.
#'
#' @param obj An object inheriting from the `AbstractKernel` class.
#' @return A string representation of the kernel.
#' @export
setGeneric("pretty_print", function(obj) standardGeneric("pretty_print"))

#' @title Generic Function for Getting Hyperparameters
#' @description
#' This generic function retrieves the hyperparameters of a kernel object.
#'
#' @param obj An object inheriting from the `AbstractKernel` class.
#' @return A list of hyperparameters.
#' @export
setGeneric("gt_HPs", function(obj) standardGeneric("gt_HPs"))

#' @title Generic Function for Kernel Derivative Computation
#' @description
#' This generic function computes the derivative of the kernel function with respect to a parameter.
#'
#' @param obj An object inheriting from the `AbstractKernel` class.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setGeneric("kernel_deriv", function(obj, x, y, param) standardGeneric("kernel_deriv"))

#' @title Generic Function for Kernel Derivative Expectation Computation
#' @description
#' This generic function computes the expected derivative of the kernel function with respect to a parameter.
#'
#' @param obj An object inheriting from the `AbstractKernel` class.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A matrix representing the expected derivative of the kernel evaluations.
#' @export
setGeneric("kernel_deriv_exp", function(obj, x, y, param) standardGeneric("kernel_deriv_exp"))

#' @title Pairwise Kernel Method for AbstractKernel
#' @description
#' This method is intended to compute the pairwise kernel matrix for an `AbstractKernel` object.
#' It must be implemented in subclasses.
#'
#' @param obj An object of class `AbstractKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @export
setMethod(
  "pairwise_kernel", "AbstractKernel",
  function(obj, x, y) {
    stop("pairwise_kernel method must be implemented in subclass")
  }
)

#' @title Show Method for AbstractKernel
#' @description
#' This method provides a basic display for objects of class `AbstractKernel`.
#'
#' @param object An object of class `AbstractKernel`.
#' @export
setMethod("show", "AbstractKernel", function(object) {
  cat("Abstract Kernel Object\n")
})

#' @title Pretty Print Method for AbstractKernel
#' @description
#' This method provides a string representation for objects of class `AbstractKernel`.
#'
#' @param obj An object of class `AbstractKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "AbstractKernel", function(obj) {
  "AbstractKernel()"
})

#' @title Hyperparameters Method for AbstractKernel
#' @description
#' This method is intended to retrieve hyperparameters from an `AbstractKernel` object.
#' It must be implemented in subclasses.
#'
#' @param obj An object of class `AbstractKernel`.
#' @export
setMethod(
  "gt_HPs", "AbstractKernel",
  function(obj) {
    stop("gt_HPs method must be implemented in subclass")
  }
)

#' @title Kernel Derivative Method for AbstractKernel
#' @description
#' This method computes the derivative of the kernel function with respect to a parameter.
#' It must be implemented in subclasses.
#'
#' @param obj An object of class `AbstractKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @export
setMethod(
  "kernel_deriv", "AbstractKernel",
  function(obj, x, y, param) {
    stop("kernel_deriv method must be implemented in subclass")
  }
)



# SumKernel ---------------------------------------------------------------
#' @title Sum Kernel Class
#' @description
#' A class representing a sum of kernels, which inherits from the `AbstractKernel` class.
#' This kernel allows combining multiple kernel objects by summing their outputs.
#'
#' @slot kernels A list of kernel objects that are summed together.
#' @export
setClass("SumKernel",
         contains = "AbstractKernel",
         slots = c(kernels = "list")
)

#' @title Initialize Method for SumKernel
#' @description
#' Initializes an instance of the `SumKernel` class with a list of kernel objects.
#'
#' @param .Object An object of class `SumKernel`.
#' @param kernels A list of kernel objects to be summed.
#' @return An initialized object of class `SumKernel`.
#' @export
setMethod(
  "initialize", "SumKernel",
  function(.Object, kernels) {
    .Object@kernels <- lapply(kernels, ensure_abstract_kernel)
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for SumKernel
#' @description
#' Computes the pairwise kernel matrix for a `SumKernel` object by summing the pairwise kernel matrices of its constituent kernels.
#'
#' @param obj An object of class `SumKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the sum of the pairwise kernel evaluations of the constituent kernels.
#' @export
setMethod(
  "pairwise_kernel", "SumKernel",
  function(obj, x, y) {
    result <- 0
    for (kernel in obj@kernels) {
      result <- result + pairwise_kernel(kernel, x, y)
    }
    return(result)
  }
)

#' @title Show Method for SumKernel
#' @description
#' Provides a display for objects of class `SumKernel`.
#'
#' @param object An object of class `SumKernel`.
#' @export
setMethod("show", "SumKernel", function(object) {
  cat("Sum Kernel:\n")
  for (i in seq_along(object@kernels)) {
    cat("  Kernel", i, ":\n")
    show(object@kernels[[i]])
  }
})

#' @title Pretty Print Method for SumKernel
#' @description
#' Provides a string representation for objects of class `SumKernel`.
#'
#' @param obj An object of class `SumKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "SumKernel", function(obj) {
  kernel_strings <- sapply(obj@kernels, pretty_print)
  paste0("[", paste(kernel_strings, collapse = " + "), "]")
})

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 ANY.
#' @return An object of class SumKernel.
#' @export
setMethod(
  "+", signature(e1 = "AbstractKernel", e2 = "ANY"),
  function(e1, e2) {
    new("SumKernel", kernels = list(e1, ensure_abstract_kernel(e2)))
  }
)

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 ANY.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class SumKernel.
#' @export
setMethod(
  "+", signature(e1 = "ANY", e2 = "AbstractKernel"),
  function(e1, e2) {
    new("SumKernel", kernels = list(ensure_abstract_kernel(e1), e2))
  }
)

#' @title Addition Method for AbstractKernel
#' @description This method defines the addition operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class SumKernel.
#' @export
setMethod(
  "+", signature(e1 = "AbstractKernel", e2 = "AbstractKernel"),
  function(e1, e2) {
    new("SumKernel", kernels = list(e1, e2))
  }
)

#' @title Hyperparameters Method for SumKernel
#' @description
#' Retrieves the hyperparameters of a `SumKernel` object by combining the hyperparameters of its constituent kernels.
#'
#' @param obj An object of class `SumKernel`.
#' @return A list of hyperparameters from all constituent kernels.
#' @export
setMethod(
  "gt_HPs", "SumKernel",
  function(obj) {
    hps <- list()
    for (kernel in obj@kernels) {
      hps <- c(hps, gt_HPs(kernel))
    }
    return(hps)
  }
)

#' @title Kernel Derivative Method for SumKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `SumKernel` object.
#' The derivative is computed as the sum of the derivatives of the constituent kernels.
#'
#' @param obj An object of class `SumKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "SumKernel",
  function(obj, x, y, param) {
    deriv <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (kernel in obj@kernels) {
      if (param %in% names(gt_HPs(kernel))) {
        deriv <- deriv + kernel_deriv(kernel, x, y, param)
      }
    }
    return(deriv)
  }
)

# ProductKernel -----------------------------------------------------------
#' @title Product Kernel Class
#' @description
#' A class representing a product of kernels, which inherits from the `AbstractKernel` class.
#' This kernel allows combining multiple kernel objects by multiplying their outputs.
#'
#' @slot kernels A list of kernel objects that are multiplied together.
#' @export
setClass("ProductKernel",
         contains = "AbstractKernel",
         slots = c(kernels = "list")
)

#' @title Initialize Method for ProductKernel
#' @description
#' Initializes an instance of the `ProductKernel` class with a list of kernel objects.
#'
#' @param .Object An object of class `ProductKernel`.
#' @param kernels A list of kernel objects to be multiplied.
#' @return An initialized object of class `ProductKernel`.
#' @export
setMethod(
  "initialize", "ProductKernel",
  function(.Object, kernels) {
    .Object@kernels <- lapply(kernels, ensure_abstract_kernel)
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for ProductKernel
#' @description
#' Computes the pairwise kernel matrix for a `ProductKernel` object by multiplying the pairwise kernel matrices of its constituent kernels.
#'
#' @param obj An object of class `ProductKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the product of the pairwise kernel evaluations of the constituent kernels.
#' @export
setMethod(
  "pairwise_kernel", "ProductKernel",
  function(obj, x, y) {
    result <- 1
    for (kernel in obj@kernels) {
      result <- result * pairwise_kernel(kernel, x, y)
    }
    return(result)
  }
)

#' @title Show Method for ProductKernel
#' @description
#' Provides a display for objects of class `ProductKernel`.
#'
#' @param object An object of class `ProductKernel`.
#' @export
setMethod("show", "ProductKernel", function(object) {
  cat("Product Kernel:\n")
  for (i in seq_along(object@kernels)) {
    cat("  Kernel", i, ":\n")
    show(object@kernels[[i]])
  }
})

#' @title Pretty Print Method for ProductKernel
#' @description
#' Provides a string representation for objects of class `ProductKernel`.
#'
#' @param obj An object of class `ProductKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "ProductKernel", function(obj) {
  kernel_strings <- sapply(obj@kernels, pretty_print)
  paste0("[", paste(kernel_strings, collapse = " * "), "]")
})

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 ANY.
#' @return An object of class ProductKernel.
#' @export
setMethod(
  "*", signature(e1 = "AbstractKernel", e2 = "ANY"),
  function(e1, e2) {
    new("ProductKernel", kernels = list(e1, ensure_abstract_kernel(e2)))
  }
)

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 ANY.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class ProductKernel.
#' @export
setMethod(
  "*", signature(e1 = "ANY", e2 = "AbstractKernel"),
  function(e1, e2) {
    new("ProductKernel", kernels = list(ensure_abstract_kernel(e1), e2))
  }
)

#' @title Multiplication Method for AbstractKernel
#' @description This method defines the multiplication operation for AbstractKernel objects.
#' @param e1 An object of class AbstractKernel.
#' @param e2 An object of class AbstractKernel.
#' @return An object of class ProductKernel.
#' @export
setMethod(
  "*", signature(e1 = "AbstractKernel", e2 = "AbstractKernel"),
  function(e1, e2) {
    new("ProductKernel", kernels = list(e1, e2))
  }
)

#' @title Hyperparameters Method for ProductKernel
#' @description
#' Retrieves the hyperparameters of a `ProductKernel` object by combining the hyperparameters of its constituent kernels.
#'
#' @param obj An object of class `ProductKernel`.
#' @return A list of hyperparameters from all constituent kernels.
#' @export
setMethod(
  "gt_HPs", "ProductKernel",
  function(obj) {
    hps <- list()
    for (kernel in obj@kernels) {
      hps <- c(hps, gt_HPs(kernel))
    }
    return(hps)
  }
)

#' @title Kernel Derivative Method for ProductKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `ProductKernel` object.
#' The derivative is computed using the product rule of differentiation.
#'
#' @param obj An object of class `ProductKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "ProductKernel",
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
  }
)


# ConstantKernel ----------------------------------------------------------
#' @title Constant Kernel Class
#' @description
#' A class representing a constant kernel, which inherits from the `AbstractKernel` class.
#' This kernel returns a constant value for all pairwise evaluations.
#'
#' @slot value_c A numeric value representing the constant value of the kernel.
#' @export
setClass("ConstantKernel",
         contains = "AbstractKernel",
         slots = c(value_c = "numeric")
)

#' @title Initialize Method for ConstantKernel
#' @description
#' Initializes an instance of the `ConstantKernel` class with a specified constant value.
#'
#' @param .Object An object of class `ConstantKernel`.
#' @param value_c A numeric value for the constant kernel. Defaults to a random uniform value between 0 and 3.
#' @return An initialized object of class `ConstantKernel`.
#' @export
setMethod(
  "initialize", "ConstantKernel",
  function(.Object, value_c = 1) {
    .Object@value_c <- value_c
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for ConstantKernel
#' @description
#' Computes the pairwise kernel matrix for a `ConstantKernel` object.
#'
#' Let \eqn{c \in \mathbb{R}_{+}^{*}}. Consider the constant kernel defined by:
#'
#' \eqn{K_{\text{cst}}(x, x') = c \cdot \begin{bmatrix}
#' 1 & 1 & \cdots & 1 \\
#' 1 & 1 & \cdots & 1 \\
#' \vdots & \vdots & \ddots & \vdots \\
#' 1 & 1 & \cdots & 1
#' \end{bmatrix}}
#'
#' The kernel matrix is filled with the constant value of the kernel.
#'
#' @param obj An object of class `ConstantKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix filled with the constant value of the kernel.
#' @export
setMethod(
  "pairwise_kernel", "ConstantKernel",
  function(obj, x, y) {
    matrix(obj@value_c, nrow = nrow(x), ncol = nrow(y))
  }
)

#' @title Show Method for ConstantKernel
#' @description
#' Provides a display for objects of class `ConstantKernel`.
#'
#' @param object An object of class `ConstantKernel`.
#' @export
setMethod("show", "ConstantKernel", function(object) {
  cat("Constant Kernel:\n")
  cat("  Value:", object@value_c, "\n")
})

#' @title Pretty Print Method for ConstantKernel
#' @description
#' Provides a string representation for objects of class `ConstantKernel`.
#'
#' @param obj An object of class `ConstantKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "ConstantKernel", function(obj) {
  sprintf("ConstantKernel(%.2f)", obj@value_c)
})

#' @title Hyperparameters Method for ConstantKernel
#' @description
#' Retrieves the hyperparameters of a `ConstantKernel` object.
#'
#' @param obj An object of class `ConstantKernel`.
#' @return A list containing the constant value of the kernel.
#' @export
setMethod(
  "gt_HPs", "ConstantKernel",
  function(obj) {
    list(value_c = obj@value_c)
  }
)

#' @title Kernel Derivative Method for ConstantKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `ConstantKernel` object.
#' The derivative with respect to the constant value is a matrix of ones.
#'
#' The derivative of \eqn{K_{\text{cst}}(x, x')} with respect to \eqn{c} is given by:
#'
#' \eqn{\frac{dK_{\text{cst}}(x, x')}{dc} = \begin{bmatrix}
#' 1 & 1 & \cdots & 1 \\
#' 1 & 1 & \cdots & 1 \\
#' \vdots & \vdots & \ddots & \vdots \\
#' 1 & 1 & \cdots & 1
#' \end{bmatrix}}
#'
#' @param obj An object of class `ConstantKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "ConstantKernel",
  function(obj, x, y, param) {
    if (param == "value_c") {
      return(matrix(1, nrow = nrow(x), ncol = nrow(y)))
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)


# NoiseKernel ------------------------------------------------------------
#' @title Noise Kernel Class
#' @description
#' A class representing a noise kernel, which inherits from the `AbstractKernel` class.
#' This kernel models noise in the data by returning a diagonal matrix scaled by a constant value.
#'
#' @slot value_c A numeric value representing the noise variance.
#' @export
setClass("NoiseKernel",
         contains = "AbstractKernel",
         slots = c(value_c = "numeric")
)

#' @title Initialize Method for NoiseKernel
#' @description
#' Initializes an instance of the `NoiseKernel` class with a specified noise variance value.
#'
#' @param .Object An object of class `NoiseKernel`.
#' @param value_c A numeric value for the noise variance. Defaults to a random uniform value between 0 and 3.
#' @return An initialized object of class `NoiseKernel`.
#' @export
setMethod(
  "initialize", "NoiseKernel",
  function(.Object, value_c = 1) {
    .Object@value_c <- value_c
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for NoiseKernel
#' @description
#' Computes the pairwise kernel matrix for a `NoiseKernel` object.
#' Let \eqn{c \in \mathbb{R}_{+}^{*}}. Consider the noise kernel defined by:
#'
#' \eqn{K_{\text{noise}}(x, x') = c \cdot I_n}
#'
#' The kernel matrix is a diagonal matrix scaled by the noise variance value.
#'
#' @param obj An object of class `NoiseKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A diagonal matrix scaled by the noise variance value.
#' @export
setMethod(
  "pairwise_kernel", "NoiseKernel",
  function(obj, x, y) {
    obj@value_c * diag(nrow(x))
  }
)

#' @title Show Method for NoiseKernel
#' @description
#' Provides a display for objects of class `NoiseKernel`.
#'
#' @param object An object of class `NoiseKernel`.
#' @export
setMethod("show", "NoiseKernel", function(object) {
  cat("Noise Kernel:\n")
  cat("  Value:", object@value_c, "\n")
})

#' @title Pretty Print Method for NoiseKernel
#' @description
#' Provides a string representation for objects of class `NoiseKernel`.
#'
#' @param obj An object of class `NoiseKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "NoiseKernel", function(obj) {
  sprintf("NoiseKernel(%.2f)", obj@value_c)
})

#' @title Hyperparameters Method for NoiseKernel
#' @description
#' Retrieves the hyperparameters of a `NoiseKernel` object.
#'
#' @param obj An object of class `NoiseKernel`.
#' @return A list containing the noise variance value.
#' @export
setMethod(
  "gt_HPs", "NoiseKernel",
  function(obj) {
    list(value_c = obj@value_c)
  }
)

#' @title Kernel Derivative Method for NoiseKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `NoiseKernel` object.
#' The derivative with respect to the noise variance is a diagonal matrix of ones.
#'
#' The derivative of \eqn{K_{\text{noise}}(x, x')} with respect to \eqn{c} is given by:
#'
#' \eqn{\frac{dK_{\text{noise}}(x, x')}{dc} = I_n = \begin{bmatrix}
#' 1 & 0 & \cdots & 0 \\
#' 0 & 1 & \cdots & 0 \\
#' \vdots & \vdots & \ddots & \vdots \\
#' 0 & 0 & \cdots & 1
#' \end{bmatrix}}
#'
#' @param obj An object of class `NoiseKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed.
#' @return A diagonal matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "NoiseKernel",
  function(obj, x, y, param) {
    if (param == "value_c") {
      return(diag(nrow(x)))
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)


# SEKernel ----------------------------------------------------------------
#' @title Squared Exponential Kernel Class
#' @description
#' A class representing a Squared Exponential (SE) kernel, which inherits from the `AbstractKernel` class.
#' The SE kernel is commonly used in Gaussian processes and is defined by a variance and length scale parameter.
#'
#' @slot variance_se A numeric value representing the signal variance.
#' @slot length_scale_se A numeric value representing the length scale parameter.
#' @export
setClass("SEKernel",
         contains = "AbstractKernel",
         slots = c(variance_se = "numeric", length_scale_se = "numeric")
)

#' @title Initialize Method for SEKernel
#' @description
#' Initializes an instance of the `SEKernel` class with specified variance and length scale parameters.
#'
#' @param .Object An object of class `SEKernel`.
#' @param variance_se A numeric value for the signal variance. Defaults to a random uniform value between 0 and 3.
#' @param length_scale_se A numeric value for the length scale parameter. Defaults to a random uniform value between 0 and 3.
#' @return An initialized object of class `SEKernel`.
#' @export
setMethod(
  "initialize", "SEKernel",
  function(.Object, variance_se = 1, length_scale_se = 1) {
    .Object@variance_se <- variance_se
    .Object@length_scale_se <- length_scale_se
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for SEKernel
#' @description
#' Computes the pairwise kernel matrix for an `SEKernel` object using the squared exponential kernel formula.
#' Let \eqn{\sigma \in \mathbb{R}_{+}} and \eqn{\ell \in \mathbb{R}_{+}^{*}}.
#' The squared exponential kernel is defined by:
#'
#' \eqn{K_{\text{SE}}(x, x') = \sigma^2 \exp\left(-\frac{\|x - x'\|^2}{2\ell^2}\right)},
#'
#' where \eqn{\sigma} is the signal standard deviation and \eqn{\ell} is the length scale parameter.
#'
#' @param obj An object of class `SEKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "SEKernel",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    return(obj@variance_se * exp(-dx / (2 * obj@length_scale_se^2)))
  }
)



#' @title Kernel Derivative Method for SEKernel
#' @description Computes the derivative of the kernel function with respect to a parameter for an `SEKernel` object.
#' Supports derivatives with respect to the variance and length scale parameters.
#'
#' The derivative of \eqn{K_{\text{SE}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{SE}}(x, x')}{d\sigma^2} = \exp\left(-\frac{\|x - x'\|^2}{2\ell^2}\right)}
#'
#' \eqn{\frac{dK_{\text{SE}}(x, x')}{d\sigma^2} = \frac{K_{\text{SE}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{SE}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \eqn{\frac{dK_{\text{SE}}(x, x')}{d\ell} = \sigma^2 \exp\left(-\frac{\|x - x'\|^2}{2\ell^2}\right) \cdot \frac{\|x - x'\|^2}{\ell^3}}
#'
#' \eqn{\frac{dK_{\text{SE}}(x, x')}{d\ell} = K_{\text{SE}}(x, x') \cdot \frac{\|x - x'\|^2}{\ell^3}}
#'
#' @param obj An object of class `SEKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Can be "variance_se" or "length_scale_se".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "SEKernel",
  function(obj, x, y, param) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    if (param == "variance_se") {
      return(pairwise_kernel(obj, x, y) / obj@variance_se)
    } else if (param == "length_scale_se") {
      return(obj@variance_se * exp(-dx / (2 * obj@length_scale_se^2)) * dx / (obj@length_scale_se^3))
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)

#' @title Show Method for SEKernel
#' @description
#' Provides a display for objects of class `SEKernel`.
#'
#' @param object An object of class `SEKernel`.
#' @export
setMethod("show", "SEKernel", function(object) {
  cat("Squared Exponential Kernel:\n")
  cat("  Variance:", object@variance_se, "\n")
  cat("  Length Scale:", object@length_scale_se, "\n")
})

#' @title Pretty Print Method for SEKernel
#' @description
#' Provides a string representation for objects of class `SEKernel`.
#'
#' @param obj An object of class `SEKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "SEKernel", function(obj) {
  sprintf("SEKernel(variance=%.2f, length_scale=%.2f)", obj@variance_se, obj@length_scale_se)
})

#' @title Hyperparameters Method for SEKernel
#' @description
#' Retrieves the hyperparameters of an `SEKernel` object.
#'
#' @param obj An object of class `SEKernel`.
#' @return A list containing the variance and length scale parameters.
#' @export
setMethod(
  "gt_HPs", "SEKernel",
  function(obj) {
    list(variance_se = obj@variance_se, length_scale_se = obj@length_scale_se)
  }
)


# LinearKernel -------------------------------------------------------------
#' @title Linear Kernel Class
#' @description
#' A class representing a linear kernel, which inherits from the `AbstractKernel` class.
#' The linear kernel is defined by two variance parameters and a constant term.
#'
#' @slot sigma2_b A numeric value representing the bias variance.
#' @slot sigma2_v A numeric value representing the variance scaling factor.
#' @slot c A numeric value representing the constant offset.
#' @export
setClass("LinearKernel",
         contains = "AbstractKernel",
         slots = c(sigma2_b = "numeric", sigma2_v = "numeric", c = "numeric")
)

#' @title Initialize Method for LinearKernel
#' @description
#' Initializes an instance of the `LinearKernel` class with specified bias variance, variance scaling factor, and constant offset.
#'
#' @param .Object An object of class `LinearKernel`.
#' @param sigma2_b A numeric value for the bias variance. Defaults to a random uniform value between 0 and 3.
#' @param sigma2_v A numeric value for the variance scaling factor. Defaults to a random uniform value between 0 and 3.
#' @param c A numeric value for the constant offset. Defaults to a random uniform value between 0 and 3.
#' @return An initialized object of class `LinearKernel`.
#' @export
setMethod(
  "initialize", "LinearKernel",
  function(.Object, sigma2_b = 1, sigma2_v = 1, c = 1) {
    .Object@sigma2_b <- sigma2_b
    .Object@sigma2_v <- sigma2_v
    .Object@c <- c
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for LinearKernel
#' @description
#' Computes the pairwise kernel matrix for a `LinearKernel` object using the linear kernel formula.
#' Let \eqn{\sigma^2_b \in \mathbb{R}_{+}}, \eqn{\sigma^2_v \in \mathbb{R}_{+}}, and \eqn{c \in \mathbb{R}}.
#' Consider the linear kernel defined by:
#'
#' \eqn{K_{\text{Lin}}(x, x') = \sigma^2_b + \sigma^2_v (x - c)(x' - c)}
#'
#' @param obj An object of class `LinearKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "LinearKernel",
  function(obj, x, y) {
    x_centered <- x - obj@c
    y_centered <- y - obj@c
    product <- tcrossprod(x_centered, y_centered)
    return(obj@sigma2_b + obj@sigma2_v * product)
  }
)

#' @title Kernel Derivative Method for LinearKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `LinearKernel` object.
#' Supports derivatives with respect to the bias variance, variance scaling factor, and constant offset.
#'
#' The derivative of \eqn{K_{\text{Lin}}(x, x')} with respect to \eqn{\sigma^2_b} is given by:
#'
#' \eqn{\frac{dK_{\text{Lin}}(x, x')}{d\sigma^2_b} = 1}
#'
#' The derivative of \eqn{K_{\text{Lin}}(x, x')} with respect to \eqn{\sigma^2_v} is given by:
#'
#' \eqn{\frac{dK_{\text{Lin}}(x, x')}{d\sigma^2_v} = (x - c)(x' - c)}
#'
#' \eqn{\frac{dK_{\text{Lin}}(x, x')}{d\sigma^2_v} = \frac{K_{\text{Lin}}(x, x')}{\sigma^2_v}}
#'
#' The derivative of \eqn{K_{\text{Lin}}(x, x')} with respect to \eqn{c} is given by:
#'
#' \eqn{\frac{dK_{\text{Lin}}(x, x')}{dc} = \sigma^2_v \left( 2c - x - x' \right)}
#'
#' \eqn{\frac{dK_{\text{Lin}}(x, x')}{dc} = -\sigma^2_v \left( (x - c) + (x' - c) \right)}
#'
#' @param obj An object of class `LinearKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Can be "sigma2_b", "sigma2_v", or "c".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "LinearKernel",
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
  }
)

#' @title Show Method for LinearKernel
#' @description
#' Provides a display for objects of class `LinearKernel`.
#'
#' @param object An object of class `LinearKernel`.
#' @export
setMethod("show", "LinearKernel", function(object) {
  cat("Linear Kernel:\n")
  cat("  Sigma squared b:", object@sigma2_b, "\n")
  cat("  Sigma squared v:", object@sigma2_v, "\n")
  cat("  c:", object@c, "\n")
})

#' @title Pretty Print Method for LinearKernel
#' @description
#' Provides a string representation for objects of class `LinearKernel`.
#'
#' @param obj An object of class `LinearKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "LinearKernel", function(obj) {
  sprintf("LinearKernel(sigma2_b=%.2f, sigma2_v=%.2f, c=%.2f)", obj@sigma2_b, obj@sigma2_v, obj@c)
})

#' @title Hyperparameters Method for LinearKernel
#' @description
#' Retrieves the hyperparameters of a `LinearKernel` object.
#'
#' @param obj An object of class `LinearKernel`.
#' @return A list containing the bias variance, variance scaling factor, and constant offset.
#' @export
setMethod(
  "gt_HPs", "LinearKernel",
  function(obj) {
    list(sigma2_b = obj@sigma2_b, sigma2_v = obj@sigma2_v, c = obj@c)
  }
)


# RationalQuadraticKernel -------------------------------------------------
#' @title Rational Quadratic Kernel Class
#' @description
#' A class representing a Rational Quadratic kernel, which inherits from the `AbstractKernel` class.
#' The Rational Quadratic kernel is a scale mixture of Squared Exponential kernels with different characteristic length-scales.
#'
#' @slot variance_rq A numeric value representing the signal variance.
#' @slot length_scale_rq A numeric value representing the length scale parameter.
#' @slot alpha_rq A numeric value representing the scale mixture parameter.
#' @export
setClass("RationalQuadraticKernel",
         contains = "AbstractKernel",
         slots = c(variance_rq = "numeric", length_scale_rq = "numeric", alpha_rq = "numeric")
)

#' @title Initialize Method for RationalQuadraticKernel
#' @description
#' Initializes an instance of the `RationalQuadraticKernel` class with specified variance, length scale, and alpha parameters.
#'
#' @param .Object An object of class `RationalQuadraticKernel`.
#' @param variance_rq A numeric value for the signal variance. Defaults to a random uniform value between 0 and 3.
#' @param length_scale_rq A numeric value for the length scale parameter. Defaults to a random uniform value between 0 and 3.
#' @param alpha_rq A numeric value for the scale mixture parameter. Defaults to a random uniform value between 0 and 3.
#' @return An initialized object of class `RationalQuadraticKernel`.
#' @export
setMethod(
  "initialize", "RationalQuadraticKernel",
  function(.Object, variance_rq = 1, length_scale_rq = 1, alpha_rq = 0.1) {
    .Object@variance_rq <- variance_rq
    .Object@length_scale_rq <- length_scale_rq
    .Object@alpha_rq <- alpha_rq
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for RationalQuadraticKernel
#' @description
#' Computes the pairwise kernel matrix for a `RationalQuadraticKernel` object using the Rational Quadratic kernel formula.
#' Let \eqn{\sigma^2 \in \mathbb{R}_{+}}, \eqn{\alpha \in \mathbb{R}_{+}^{*}}, and \eqn{\ell \in \mathbb{R}_{+}^{*}}.
#' Consider the RQ kernel defined by:
#'
#' \eqn{K_{\text{RQ}}(x, x') = \sigma^2 \left(1 + \frac{\|x - x'\|^2}{2 \alpha \ell^2}\right)^{-\alpha}}
#'
#' @param obj An object of class `RationalQuadraticKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "RationalQuadraticKernel",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    return(obj@variance_rq * (1 + dx / (2 * obj@alpha_rq * obj@length_scale_rq^2))^(-obj@alpha_rq))
  }
)

#' @title Kernel Derivative Method for RationalQuadraticKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `RationalQuadraticKernel` object.
#' Supports derivatives with respect to the variance, length scale, and alpha parameters.
#'
#' The derivative of \eqn{K_{\text{RQ}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\sigma^2} = \left(1 + \frac{\|x - x'\|^2}{2 \alpha \ell^2}\right)^{-\alpha}}
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\sigma^2} = \frac{K_{\text{RQ}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{RQ}}(x, x')} with respect to \eqn{\alpha} is given by:
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\alpha} = \sigma^2 \left( \frac{\|x - x'\|^2}{2 \alpha \ell^2} + 1 \right)^{-\alpha} \left( \frac{\|x - x'\|^2}{2\alpha \ell^2} \left( \frac{\|x - x'\|^2}{2 \alpha \ell^2} + 1 \right)^{-1} - \ln \left( \frac{\|x - x'\|^2}{2\alpha \ell^2} + 1 \right) \right)}
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\alpha} = -\sigma^2 \left( \frac{\|x - x'\|^2}{2 \alpha \ell^2} + 1 \right)^{-\alpha} \left( \frac{(2 \alpha \ell^2 + \|x - x'\|^2) \ln \left( \frac{\|x - x'\|^2}{2 \alpha \ell^2} + 1 \right) - \|x - x'\|^2}{2 \alpha \ell^2 + \|x - x'\|^2} \right)}
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\alpha} = -K_{\text{RQ}}(x, x') \left( \frac{(2 \alpha \ell^2 + \|x - x'\|^2) \ln \left( \frac{\|x - x'\|^2}{2 \alpha \ell^2} + 1 \right) - \|x - x'\|^2}{2 \alpha \ell^2 + \|x - x'\|^2} \right)}
#'
#' The derivative of \eqn{K_{\text{RQ}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\ell} = \sigma^2 \cdot \left(1 + \frac{\|x - x'\|^2}{2 \alpha \ell^2}\right)^{-\alpha-1} \cdot \left( \frac{\|x - x'\|^2}{\ell^3} \right)}
#'
#' \eqn{\frac{dK_{\text{RQ}}(x, x')}{d\ell} = K_{\text{RQ}}(x, x') \cdot \frac{\|x - x'\|^2}{\ell^3 \left(1 + \frac{\|x - x'\|^2}{2 \alpha \ell^2}\right)}}
#'
#' @param obj An object of class `RationalQuadraticKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Can be "variance_rq", "length_scale_rq", or "alpha_rq".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "RationalQuadraticKernel",
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
  }
)

#' @title Show Method for RationalQuadraticKernel
#' @description
#' Provides a display for objects of class `RationalQuadraticKernel`.
#'
#' @param object An object of class `RationalQuadraticKernel`.
#' @export
setMethod("show", "RationalQuadraticKernel", function(object) {
  cat("Rational Quadratic Kernel:\n")
  cat("  Variance:", object@variance_rq, "\n")
  cat("  Length Scale:", object@length_scale_rq, "\n")
  cat("  Alpha:", object@alpha_rq, "\n")
})

#' @title Pretty Print Method for RationalQuadraticKernel
#' @description
#' Provides a string representation for objects of class `RationalQuadraticKernel`.
#'
#' @param obj An object of class `RationalQuadraticKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "RationalQuadraticKernel", function(obj) {
  sprintf("RationalQuadraticKernel(variance=%.2f, length_scale=%.2f, alpha=%.2f)", obj@variance_rq, obj@length_scale_rq, obj@alpha_rq)
})

#' @title Hyperparameters Method for RationalQuadraticKernel
#' @description
#' Retrieves the hyperparameters of a `RationalQuadraticKernel` object.
#'
#' @param obj An object of class `RationalQuadraticKernel`.
#' @return A list containing the variance, length scale, and alpha parameters.
#' @export
setMethod(
  "gt_HPs", "RationalQuadraticKernel",
  function(obj) {
    list(variance_rq = obj@variance_rq, length_scale_rq = obj@length_scale_rq, alpha_rq = obj@alpha_rq)
  }
)


# PeriodicKernel ----------------------------------------------------------
#' @title Periodic Kernel Class
#' @description
#' A class representing a Periodic kernel, which inherits from the `AbstractKernel` class.
#' The Periodic kernel is used to model periodic functions and is defined by a variance, length scale, and period.
#'
#' @slot variance_per A numeric value representing the signal variance.
#' @slot length_scale_per A numeric value representing the length scale parameter.
#' @slot period A numeric value representing the period of the periodic function.
#' @export
setClass("PeriodicKernel",
         contains = "AbstractKernel",
         slots = c(variance_per = "numeric", length_scale_per = "numeric", period = "numeric")
)

#' @title Initialize Method for PeriodicKernel
#' @description
#' Initializes an instance of the `PeriodicKernel` class with specified variance, length scale, and period parameters.
#'
#' @param .Object An object of class `PeriodicKernel`.
#' @param variance_per A numeric value for the signal variance.
#' @param length_scale_per A numeric value for the length scale parameter.
#' @param period A numeric value for the period.
#' @return An initialized object of class `PeriodicKernel`.
#' @export
setMethod(
  "initialize", "PeriodicKernel",
  function(.Object, variance_per =1, length_scale_per = 1, period = 1) {
    .Object@variance_per <- variance_per
    .Object@length_scale_per <- length_scale_per
    .Object@period <- period
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for PeriodicKernel
#' @description
#' Computes the pairwise kernel matrix for a `PeriodicKernel` object using the periodic kernel formula.
#' Let \eqn{\sigma^2 \in \mathbb{R}_{+}}, \eqn{\ell \in \mathbb{R}}, and \eqn{p \in \mathbb{R}}.
#' Consider the periodic kernel defined by:
#'
#' \eqn{K_{\text{Per}}(x, x') = \sigma^2 \exp \left( -\frac{2}{\ell^2} \sin^2 \left( \pi \frac{x - x'}{p} \right) \right)}
#'
#' @param obj An object of class `PeriodicKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "PeriodicKernel",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    return(obj@variance_per * exp(-2 * (sin(pi * dx / obj@period)^2) / obj@length_scale_per^2))
  }
)

#' @title Kernel Derivative Method for PeriodicKernel
#' @description
#' Computes the derivative of the kernel function with respect to a parameter for a `PeriodicKernel` object.
#' Supports derivatives with respect to the variance, length scale, and period parameters.
#'
#' The derivative of \eqn{K_{\text{Per}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{Per}}(x, x')}{d\sigma^2} = \exp \left( -\frac{2}{\ell^2} \sin^2 \left( \pi \frac{x - x'}{p} \right) \right)}
#'
#' \eqn{\frac{dK_{\text{Per}}(x, x')}{d\sigma^2} = \frac{K_{\text{Per}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{Per}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \eqn{\frac{dK_{\text{Per}}(x, x')}{d\ell} = \frac{4\sigma^2 \sin^2\left(\frac{\pi(x-x')}{p}\right) \exp\left(-\frac{2\sin^2\left(\frac{\pi(x-x')}{p}\right)}{\ell^2}\right)}{\ell^3}}
#'
#' The derivative of \eqn{K_{\text{Per}}(x, x')} with respect to \eqn{p} is given by:
#'
#' \eqn{\frac{dK_{\text{Per}}(x, x')}{dp} = \frac{4\pi\sigma^2(x-x')\sin\left(\frac{\pi(x-x')}{p}\right)\cos\left(\frac{\pi(x-x')}{p}\right)\exp\left(-\frac{2\sin^2\left(\frac{\pi(x-x')}{p}\right)}{\ell^2}\right)}{\ell^2 p^2}}
#'
#' \eqn{\frac{dK_{\text{Per}}(x, x')}{dp} = \frac{2\pi\sigma^2(x-x')\sin\left(\frac{2\pi(x-x')}{p}\right)\exp\left(-\frac{2\sin^2\left(\frac{\pi(x-x')}{p}\right)}{\ell^2}\right)}{\ell^2 p^2}}
#'
#' @param obj An object of class `PeriodicKernel`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Can be "variance_per", "length_scale_per", or "period".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "PeriodicKernel",
  function(obj, x, y, param) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    if (param == "variance_per") {
      return(exp(-2 * (sin(pi * dx / obj@period)^2) / obj@length_scale_per^2) / obj@variance_per)
    } else if (param == "length_scale_per") {
      deriv <- 4 * obj@variance_per * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 * exp(-2 * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 / obj@length_scale_per^2) / obj@length_scale_per^3
      return(deriv[, 1, , drop = TRUE])
    } else if (param == "period") {
      deriv <- 2 * pi * obj@variance_per * outer(x, y, FUN = "-") * sin(2 * pi * outer(x, y, FUN = "-") / obj@period) * exp(-2 * sin(pi * outer(x, y, FUN = "-") / obj@period)^2 / obj@length_scale_per^2) / (obj@length_scale_per^2 * obj@period^2)
      return(deriv[, 1, , drop = TRUE])
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)

#' @title Show Method for PeriodicKernel
#' @description
#' Provides a display for objects of class `PeriodicKernel`.
#'
#' @param object An object of class `PeriodicKernel`.
#' @export
setMethod("show", "PeriodicKernel", function(object) {
  cat("Periodic Kernel:\n")
  cat("  Variance:", object@variance_per, "\n")
  cat("  Length Scale:", object@length_scale_per, "\n")
  cat("  Period:", object@period, "\n")
})

#' @title Pretty Print Method for PeriodicKernel
#' @description
#' Provides a string representation for objects of class `PeriodicKernel`.
#'
#' @param obj An object of class `PeriodicKernel`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "PeriodicKernel", function(obj) {
  sprintf("PeriodicKernel(variance=%.2f, length_scale=%.2f, period=%.2f)", obj@variance_per, obj@length_scale_per, obj@period)
})

#' @title Hyperparameters Method for PeriodicKernel
#' @description
#' Retrieves the hyperparameters of a `PeriodicKernel` object.
#'
#' @param obj An object of class `PeriodicKernel`.
#' @return A list containing the variance, length scale, and period parameters.
#' @export
setMethod(
  "gt_HPs", "PeriodicKernel",
  function(obj) {
    list(variance_per = obj@variance_per, length_scale_per = obj@length_scale_per, period = obj@period)
  }
)


# MaternKernel 1/2 ------------------------------------------------------------
#' @title Matern Kernel 1/2 Class
#' @description
#' A class representing a Matern kernel with smoothness parameter \eqn{\nu} = 1/2 , which inherits from the `AbstractKernel` class.
#' This kernel is a type of covariance function used in Gaussian process regression.
#'
#' @slot length_scale_mat A numeric value representing the length scale parameter.
#' @export
setClass("MaternKernel12",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric")
)

#' @title Initialize Method for MaternKernel12
#' @description
#' Initializes an instance of the `MaternKernel12` class with a specified length scale parameter.
#'
#' @param .Object An object of class `MaternKernel12`.
#' @param length_scale_mat A numeric value for the length scale parameter.
#' @return An initialized object of class `MaternKernel12`.
#' @export
setMethod(
  "initialize", "MaternKernel12",
  function(.Object, length_scale_mat = 1) {
    .Object@length_scale_mat <- length_scale_mat
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for MaternKernel12
#' @description
#' Computes the pairwise kernel matrix for a `MaternKernel12` object using the Matern kernel formula with \eqn{\nu = 1/2}.
#' Let \(\eqn{\ell \in \mathbb{R}_{+}^{*}}\) and \(\eqn{\sigma^2 \in \mathbb{R}_{+}}\). Consider the Matern Kernel with \eqn{\nu = 1/2} defined by:
#'
#' \eqn{K_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x') = \sigma^2 \cdot \exp\left(-\frac{\|x - x'\|^2}{\ell}\right)}
#'
#' @param obj An object of class `MaternKernel12`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "MaternKernel12",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    return(exp(-dx / obj@length_scale_mat))
  }
)

#' @title Kernel Derivative Method for MaternKernel12
#' @description
#' Computes the derivative of the kernel function with respect to the length scale parameter for a `MaternKernel12` object.
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}{d\sigma^2} = \exp\left(-\frac{\|x - x'\|^2}{\ell}\right)}
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}{d\sigma^2} = \frac{K_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}{d\ell} = -\exp\left(-\frac{\|x - x'\|^2}{\ell}\right) \cdot \|x - x'\|^2 \cdot (-\ell^{-2})}
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}{d\ell} = \frac{\|x - x'\|^2}{\ell^2} \cdot K_{\text{Mat}_{\nu = \frac{1}{2}}}(x, x')}
#'
#' @param obj An object of class `MaternKernel12`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Should be "length_scale_mat".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "MaternKernel12",
  function(obj, x, y, param) {
    if (param == "length_scale_mat") {
      dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
      return((dx / obj@length_scale_mat^2) * exp(-dx / obj@length_scale_mat))
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)

#' @title Show Method for MaternKernel12
#' @description
#' Provides a display for objects of class `MaternKernel12`.
#'
#' @param object An object of class `MaternKernel12`.
#' @export
setMethod("show", "MaternKernel12", function(object) {
  cat("Matern Kernel 1/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

#' @title Pretty Print Method for MaternKernel12
#' @description
#' Provides a string representation for objects of class `MaternKernel12`.
#'
#' @param obj An object of class `MaternKernel12`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "MaternKernel12", function(obj) {
  sprintf("MaternKernel12(length_scale=%.2f)", obj@length_scale_mat)
})

#' @title Hyperparameters Method for MaternKernel12
#' @description
#' Retrieves the hyperparameters of a `MaternKernel12` object.
#'
#' @param obj An object of class `MaternKernel12`.
#' @return A list containing the length scale parameter.
#' @export
setMethod(
  "gt_HPs", "MaternKernel12",
  function(obj) {
    list(length_scale_mat = obj@length_scale_mat)
  }
)


# MaternKernel 3/2 ------------------------------------------------------------
#' @title Matern Kernel 3/2 Class
#' @description
#' A class representing a Matern kernel with smoothness parameter \eqn{\nu = 3/2}, which inherits from the `AbstractKernel` class.
#' This kernel is a type of covariance function used in Gaussian process regression.
#'
#' @slot length_scale_mat A numeric value representing the length scale parameter.
#' @export
setClass("MaternKernel32",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric")
)

#' @title Initialize Method for MaternKernel32
#' @description
#' Initializes an instance of the `MaternKernel32` class with a specified length scale parameter.
#'
#' @param .Object An object of class `MaternKernel32`.
#' @param length_scale_mat A numeric value for the length scale parameter.
#' @return An initialized object of class `MaternKernel32`.
#' @export
setMethod(
  "initialize", "MaternKernel32",
  function(.Object, length_scale_mat = 1) {
    .Object@length_scale_mat <- length_scale_mat
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for MaternKernel32
#' @description
#' Computes the pairwise kernel matrix for a `MaternKernel32` object using the Matern kernel formula with \eqn{\nu = 3/2}.
#' Let \eqn{\ell \in \mathbb{R}_{+}^{*}} and \eqn{\sigma^2 \in \mathbb{R}_{+}}. Consider the Matern Kernel with \eqn{\nu = 3/2} defined by:
#'
#' \eqn{K_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x') = \sigma^2 \cdot \left( \frac{\sqrt{3} \cdot \|x - x'\|^2}{\ell} + 1 \right) \exp\left(-\frac{\sqrt{3} \cdot \|x - x'\|^2}{\ell}\right)}
#'
#' @param obj An object of class `MaternKernel32`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "MaternKernel32",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    sqrt3_r_div_l <- (sqrt(3) * dx) / obj@length_scale_mat
    return((1 + sqrt3_r_div_l) * exp(-sqrt3_r_div_l))
  }
)

#' @title Kernel Derivative Method for MaternKernel32
#' @description
#' Computes the derivative of the kernel function with respect to the length scale parameter for a `MaternKernel32` object.
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')}{d\sigma^2} = \left( \frac{\sqrt{3} \cdot \|x - x'\|^2}{\ell} + 1 \right) \exp\left(-\frac{\sqrt{3} \cdot \|x - x'\|^2}{\ell}\right)}
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')}{d\sigma^2} = \frac{K_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{3}{2}}}(x, x')}{d\ell} = \frac{3 \cdot (\|x - x'\|^2)^2 \cdot \exp\left(-\frac{\sqrt{3} \|x - x'\|^2}{\ell}\right)}{\ell^3}}
#'
#' @param obj An object of class `MaternKernel32`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Should be "length_scale_mat".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "MaternKernel32",
  function(obj, x, y, param) {
    if (param == "length_scale_mat") {
      dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
      return(3 * dx^2 * exp(-sqrt(3) * dx / obj@length_scale_mat) / obj@length_scale_mat^3)
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)

#' @title Show Method for MaternKernel32
#' @description
#' Provides a display for objects of class `MaternKernel32`.
#'
#' @param object An object of class `MaternKernel32`.
#' @export
setMethod("show", "MaternKernel32", function(object) {
  cat("Matern Kernel 3/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

#' @title Pretty Print Method for MaternKernel32
#' @description
#' Provides a string representation for objects of class `MaternKernel32`.
#'
#' @param obj An object of class `MaternKernel32`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "MaternKernel32", function(obj) {
  sprintf("MaternKernel32(length_scale=%.2f)", obj@length_scale_mat)
})

#' @title Hyperparameters Method for MaternKernel32
#' @description
#' Retrieves the hyperparameters of a `MaternKernel32` object.
#'
#' @param obj An object of class `MaternKernel32`.
#' @return A list containing the length scale parameter.
#' @export
setMethod(
  "gt_HPs", "MaternKernel32",
  function(obj) {
    list(length_scale_mat = obj@length_scale_mat)
  }
)

# MaternKernel 5/2 ------------------------------------------------------------
#' @title Matern Kernel 5/2 Class
#' @description
#' A class representing a Matern kernel with smoothness parameter \eqn{\nu = 5/2}, which inherits from the `AbstractKernel` class.
#' This kernel is a type of covariance function used in Gaussian process regression.
#'
#' @slot length_scale_mat A numeric value representing the length scale parameter.
#' @export
setClass("MaternKernel52",
         contains = "AbstractKernel",
         slots = c(length_scale_mat = "numeric")
)

#' @title Initialize Method for MaternKernel52
#' @description
#' Initializes an instance of the `MaternKernel52` class with a specified length scale parameter.
#'
#' @param .Object An object of class `MaternKernel52`.
#' @param length_scale_mat A numeric value for the length scale parameter.
#' @return An initialized object of class `MaternKernel52`.
#' @export
setMethod(
  "initialize", "MaternKernel52",
  function(.Object, length_scale_mat = 1) {
    .Object@length_scale_mat <- length_scale_mat
    return(.Object)
  }
)

#' @title Pairwise Kernel Method for MaternKernel52
#' @description
#' Computes the pairwise kernel matrix for a `MaternKernel52` object using the Matern kernel formula with \eqn{\nu = 5/2}.
#' Let \eqn{\ell \in \mathbb{R}_{+}^{*}} and \eqn{\sigma^2 \in \mathbb{R}_{+}}. Consider the Matern Kernel with \eqn{\nu = 5/2} defined by:
#'
#' \eqn{K_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x') = \sigma^2 \cdot \left( \left( \frac{5(\|x - x'\|^2)^2}{3\ell^2} + \frac{\sqrt{5}\|x - x'\|^2}{\ell} + 1 \right) \exp \left( -\frac{\sqrt{5}\|x - x'\|^2}{\ell} \right) \right)}
#'
#' @param obj An object of class `MaternKernel52`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @return A matrix representing the pairwise kernel evaluations.
#' @export
setMethod(
  "pairwise_kernel", "MaternKernel52",
  function(obj, x, y) {
    dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
    sqrt5_r_div_l <- (sqrt(5) * dx) / obj@length_scale_mat
    return((1 + sqrt5_r_div_l + (5.0 / 3.0) * (dx / obj@length_scale_mat)^2) * exp(-sqrt5_r_div_l))
  }
)

#' @title Kernel Derivative Method for MaternKernel52
#' @description
#' Computes the derivative of the kernel function with respect to the length scale parameter for a `MaternKernel52` object.
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')} with respect to \eqn{\sigma^2} is given by:
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')}{d\sigma^2} = \left( \frac{5(\|x - x'\|^2)^2}{3\ell^2} + \frac{\sqrt{5}\|x - x'\|^2}{\ell} + 1 \right) \exp\left(-\frac{\sqrt{5}\|x - x'\|^2}{\ell}\right)}
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')}{d\sigma^2} = \frac{K_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')}{\sigma^2}}
#'
#' The derivative of \eqn{K_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')} with respect to \eqn{\ell} is given by:
#'
#' \deqn{\frac{dK_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')}{d\ell} = \frac{\sqrt{5}\|x - x'\|^2 \exp\left(-\frac{\sqrt{5}\|x - x'\|^2}{\ell}\right) \left( \frac{5(\|x - x'\|^2)^2}{3\ell^2} + \frac{\sqrt{5}\|x - x'\|^2}{\ell} + 1 \right)}{\ell^2}}{+ \frac{\exp\left(-\frac{\sqrt{5}\|x - x'\|^2}{\ell}\right) \left( -\frac{10(\|x - x'\|^2)^2}{3\ell^3} - \frac{\sqrt{5}\|x - x'\|^2}{\ell^2} \right)}{\ell}}
#'
#' \eqn{\frac{dK_{\text{Mat}_{\nu = \frac{5}{2}}}(x, x')}{d\ell} = \frac{5(\|x - x'\|^2)^2 \exp\left(-\frac{\sqrt{5}\|x - x'\|^2}{\ell}\right) \left( \ell + \sqrt{5}\|x - x'\|^2 \right)^3}{\ell^4}}
#'
#' @param obj An object of class `MaternKernel52`.
#' @param x A matrix of input vectors.
#' @param y A matrix of input vectors.
#' @param param The parameter with respect to which the derivative is computed. Should be "length_scale_mat".
#' @return A matrix representing the derivative of the kernel evaluations.
#' @export
setMethod(
  "kernel_deriv", "MaternKernel52",
  function(obj, x, y, param) {
    if (param == "length_scale_mat") {
      dx <- outer(rowSums(x^2), rowSums(y^2), FUN = "+") - 2 * tcrossprod(x, y)
      return(5 * dx^2 * exp(-sqrt(5) * dx / obj@length_scale_mat) * (obj@length_scale_mat + sqrt(5) * dx)^3 / obj@length_scale_mat^4)
    } else {
      stop("Unknown parameter for derivative calculation.")
    }
  }
)

#' @title Show Method for MaternKernel52
#' @description
#' Provides a display for objects of class `MaternKernel52`.
#'
#' @param object An object of class `MaternKernel52`.
#' @export
setMethod("show", "MaternKernel52", function(object) {
  cat("Matern Kernel 5/2:\n")
  cat("  Length Scale:", object@length_scale_mat, "\n")
})

#' @title Pretty Print Method for MaternKernel52
#' @description
#' Provides a string representation for objects of class `MaternKernel52`.
#'
#' @param obj An object of class `MaternKernel52`.
#' @return A string representation of the object.
#' @export
setMethod("pretty_print", "MaternKernel52", function(obj) {
  sprintf("MaternKernel52(length_scale=%.2f)", obj@length_scale_mat)
})

#' @title Hyperparameters Method for MaternKernel52
#' @description
#' Retrieves the hyperparameters of a `MaternKernel52` object.
#'
#' @param obj An object of class `MaternKernel52`.
#' @return A list containing the length scale parameter.
#' @export
setMethod(
  "gt_HPs", "MaternKernel52",
  function(obj) {
    list(length_scale_mat = obj@length_scale_mat)
  }
)





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

  if (length(mean) == 1) {
    mean <- rep(mean, nrow(db))
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
optim_hp <- function(hp, db, mean, kern, post_cov, pen_diag=1e-6, verbose = FALSE) {
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
