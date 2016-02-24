# kernel_comp

#' @name composition
#' @rdname composition
#'
#' @title Compositional kernels
#' 
#' @description Construct a new kernel by combining existing kernels, either 
#' by summation, multiplication or the kronecker product.
#' 
#' Summation and multiplication require that the covariance matrices produced
#' the two kernels have the same dimension (same number of rows in the input
#' columns) and result in a kernel with the same dimension as its inputs.
#' 
#' The kronecker product doesn't require that the input functions have the 
#' same dimension, and the dimension of the output is the product of the 
#' dimensions of the inputs (i.e. an m-by-m matrix kroneckered with an 
#' n-by-n matrix gives rise to an nm-by-nm matrix).
#' 
#'
#' @template kco
#' 
#' @param \dots several kernel objects to be combined
#' @param kernel,kernel1,kernel2,X,Y kernel objects to be combined
#' @param na.rm an unused argument for consistency with the generic sum
#'  function
#' @examples
#' 
#' # construct a kernel with one feature
#' k1 <- rbf('temperature')
#' 
#' # and another with two features
#' k2 <- rbf(c('temperature', 'pressure'))
#' 
NULL

# underlying compositional kernel function
comp <- function (kernel1, kernel2, type) {
  
  # create a compositional kernel
  # type must be one of 'sum', 'prod' or 'kron'
  
  # create kernel data object
  object <- list(type = type,
                 kernel1 = kernel1,
                 kernel2 = kernel2)
  
  # create a function to evaluate it
  ans <- function(data, newdata = NULL, diag = FALSE) {
    
    compEval(object, data, newdata, type, diag = diag)    
    
  }
  
  # tell this function it's now a kernel
  ans <- as.kernel(ans)
  
  return(ans)
  
}

# evaluate compositional kernel
compEval <- function(object, data, newdata, operation, diag) {
  
  stopifnot(operation %in% c('sum', 'prod', 'kron'))
  
  # evaluate the sub kernels
  covmat1 <- object$kernel1(data, newdata, diag)
  covmat2 <- object$kernel2(data, newdata, diag)
  
  # get the compositional covariance matrix
  if (operation == 'sum') {
    
    covmat <- covmat1 + covmat2
    
  } else if (operation == 'prod') {
    
    covmat <- covmat1 * covmat2
    
  } else {
    
    covmat <- covmat1 %x% covmat2
    
  }
  
  # return the covariance matrix
  return (covmat)
  
}
