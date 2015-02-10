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
  
  comp(kernel1,
       kernel2,
       'sum')
  
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
  
  comp(kernel1,
       kernel2,
       'prod')
  
}

#' @rdname composition
#' @export
#' @examples
#' 
#' # get the kronecker product of two kernels
#' k_kron <- kron(k1, k2)
#' 
#' # evaluate the function and look at the matrix
#' image(k_kron(pressure))
#'  
kron <- function (kernel1,
                  kernel2) {
  
  # if both things are kernels, get a kronecker product kernel
  if (is.kernel(kernel1) &
        is.kernel(kernel2)) {
    
    ans <- comp(kernel1,
                kernel2,
                'kron')
    
  } else {
    
    stop ('kernel1 and kernel2 must be objects of class kernel')
    
  }
  
  # return the result
  return (ans)
  
}


#' @rdname composition
#' @export
#' @examples
#' 
#' # get the kronecker product of two kernels again
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