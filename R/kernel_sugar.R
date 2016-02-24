# functions to handle syntactic sugar and overload existing functions

# ~~~
# overload the arithmetic operators for addition, multiplication
# and the kronecker product

#' @rdname composition
#' @export
#' @examples
#'
#' # sum lots of kernels
#' k_sum <- sum(k1, k2, k1)
#' 
#' # evaluate the function and look at the matrix
#' image(k_sum(pressure))
#'  
sum.kernel <- function (..., na.rm = FALSE) {
  
  # parse all arguments as a list
  dots <- list(...)
  
  # check they're all kernels
  lapply(dots, checkKernel)

  # get length of the list
  narg <- length(dots)
  
  if (narg == 0) {
    
    # if it's empty, return NULL
    ans <- NULL
    
  } else if (narg == 1) {
    
    # if there's only one element, return it
    ans <- dots[[1]]
    
  } else {
    
    # otherwise, sum the first two
    ans <- dots[[1]] + dots[[2]]
    
  }
  
  if (narg > 2) {
    
    # if there are more than this, add the others sequentially
    for (i in 3:narg) {
      ans <- ans + dots[[i]]
    }
  }
  
  return (ans)
  
}

#' @rdname composition
#' @export
#' @examples
#' 
#' # sum two kernels
#' k_1p2 <- k1 + k2
#' 
#' # evaluate the function and look at the matrix
#' image(k_1p2(pressure))
#'  
`+.kernel` <- function (kernel1,
                        kernel2) {

  # check they're both kernels
  checkKernel(kernel1)
  checkKernel(kernel2)
  
  ans <- comp(kernel1,
       kernel2,
       'sum')
  
  return (ans)
  
}


#' @rdname composition
#' @export
#' @examples
#' 
#' # multiply lots of kernels
#' k_prod <- prod(k1, k2, k1)
#' 
#' # evaluate the function and look at the matrix
#' image(k_prod(pressure))
#'  
prod.kernel <- function (..., na.rm = FALSE) {
  
  # parse all arguments as a list
  dots <- list(...)

  # check they're all kernels
  lapply(dots, checkKernel)
  
  # get length of the list
  narg <- length(dots)
  
  if (narg == 0) {
    
    # if it's empty, return NULL
    ans <- NULL
    
  } else if (narg == 1) {
    
    # if there's only one element, return it
    ans <- dots[[1]]
    
  } else {
    
    # otherwise, sum the first two
    ans <- dots[[1]] * dots[[2]]
    
  }
  
  if (narg > 2) {
    
    # if there are more than this, add the others sequentially
    for (i in 3:narg) {
      ans <- ans * dots[[i]]
    }
  }
  
  return (ans)
  
}


#' @rdname composition
#' @export
#' @examples
#' 
#' # multiply two kernels
#' k_1_2 <- k1 * k2
#' 
#' # evaluate the function and look at the matrix
#' image(k_1_2(pressure))
#'  
`*.kernel` <- function (kernel1,
                        kernel2) {

  # check they're both kernels
  checkKernel(kernel1)
  checkKernel(kernel2)
  
  ans <- comp(kernel1,
       kernel2,
       'prod')
  
  return (ans)
  
}

#' @rdname composition
#' @export
#' 
#' @param power an integer (or integer-esque numeric) giving the power to which
#' to raise the kernel function. If \code{power} is not integer-esque (that is
#' \code{power != round(power)}) a warning is issued and the power rounded.
#' 
#' @examples
#' # get a cubic kernel
#' k <- int() + lin('pressure', c = 400, sigma = 0.003)
#' k_cu <- k ^ 3
#' 
#' # evaluate the function and look at the matrix
#' image(k_cu(pressure))
#' 
#' # look at example draws from the original and cubed kernel
#' demoKernel(k, pressure) 
#' demoKernel(k_cu, pressure) 
#'  
`^.kernel` <- function (kernel, power) {
  
  # check the kernel object
  checkKernel(kernel)
  
  # check the power is (essentially) an integer
  if (!power == round(power)) {
    
    # if not, coerce it into an integer
    power <- round(power)
    
    # and issue a warning
    warning (paste0('power must be an integer, but a non-integer was passed so\
                  power has been changed to ',
                    power))
    
  }
  
  # get a list of the same kernel
  kernel_list <- replicate(power, kernel, simplify = FALSE)
  
  # get the product of these
  ans <- do.call(prod, kernel_list)
  
  # return
  return (ans)
  
} 

setClass("kernel")

#' @rdname composition
#' @export
#' @param FUN the operation to use in the kronecker operation, included only
#'  for compatability, FUN = "*" is used regardless of what's specified
#' @param make.dimnames included only for compatability, ignored here
#' 
#' @examples
#' 
#' # get the kronecker product of two kernels
#' k_kron <- kronecker(k1, k2)
#' 
#' # evaluate the function and look at the matrix
#' image(k_kron(pressure))
#'  
setMethod(f = "kronecker", 
          signature(X = "kernel", Y = "kernel", FUN = "ANY", make.dimnames = "ANY"),
          definition = function(X, Y, ...) {
            comp(X,
                 Y,
                 'kron')
          }
)

#' @rdname composition
#' @export
#' @examples
#' 
#' # get the kronecker product of two kernels more elegantly
#' k_kron <- k1 %x% k2
#' 
#' # evaluate the function and look at the matrix
#' image(k_kron(pressure))
#'  
`%x%` <- function (kernel1,
                   kernel2) {
  
  if (is.kernel(kernel1) &
        is.kernel(kernel2)) {
    
    ans <- comp(kernel1,
                kernel2,
                'kron')
    
  } else {
    
    # otherwise, use the base kronecker function
    ans <- kronecker(kernel1,
                     kernel2)
    
  }
  
  return (ans)
  
} 

# define generic methods for kernel multiplication and division,
# export but don't document them

#' @export
`-.kernel` <- function (kernel1,
                        kernel2) {
  
  stop ('kernel subtraction is not possible as the resulting kernel \
        may well not be positive definite')
  
}

#' @export
`/.kernel` <- function (kernel1,
                        kernel2) {
  
  stop ('kernel division is not possible as the resulting kernel \
        may well not be positive definite')
  
}