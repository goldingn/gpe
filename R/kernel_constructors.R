# kernel constructor functions

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
#' @param \dots kernel objects to be combined
#' @param kernel1 A kernel object
#' @param kernel2 Another kernel object to combine it with
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
kernel.comp <- function (kernel1, kernel2, type) {
  
  # create a compositional kernel
  # type must be one of 'sum', 'prod' or 'kron'
  
  # create kernel data object
  object <- list(type = type,
                 kernel1 = kernel1,
                 kernel2 = kernel2)
  
  # create a function to evaluate it
  ans <- function(data, newdata = NULL) {
    
    evalKernel(object, data, newdata)    
    
  }
  
  # tell this function it's now a kernel
  ans <- as.kernel(ans)
  
  return(ans)
  
}


# ~~~~~~~~~~
# basis kernels

#' @title Radial basis function kernel
#' 
#' @description Construct a radial basis function (A.K.A. squared-exponential) kernel.
#' 
#' @details The rbf kernel takes the form:
#' \deqn{k_{rbf}(\mathbf{x}, \mathbf{x}') = \sigma^2 exp\left[-\frac{1}{2} {\sum\limits_{d=1}^D \left(\frac{(x_d - x_d')}{2l_d^2}\right)}^2\right]}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active, \eqn{l_d} 
#' are the characteristic lengthscales for each covariate (column) \eqn{x_d} 
#' and \eqn{\sigma^2} is the overall variance.
#' 
#' Larger values of \eqn{l_i} correspond to functions in which change less 
#' rapidly over the values of the covariates.
#' 
#' @template kco
#' @template kco_basis
#' @export
#' @name rbf
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- rbf('temperature')
#' 
#' # and another with two features
#' k2 <- rbf(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
rbf <- function (columns) {
  
  # construct an rbf kernel
  
  # create the model object and initialize parameters
  object <- list(type = 'rbf',
                 columns = columns,
                 parameters = list(sigma = 1,
                                   l = rep(1,
                                           length(columns))))
  
  # create a function to return
  ans <- function(data,
                  newdata = NULL) {
    
    evalKernel(object,
               data,
               newdata)    
    
  }
  
  # tell this function it's now a kernel
  ans <- as.kernel(ans)
  
  return(ans)
  
}
