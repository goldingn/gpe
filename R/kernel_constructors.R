# kernel constructor functions
# compositional kernel function - should this be exported?

# @title Compositional kernel
# 
# @description Create a new kernel object by combining two existing kernel objects, 
# either by summation, multiplication or convolution.
# 
# @template kco
# @template kco_comp
# 
# @param type A string giving the type of convolution to perform, either 
# 'sum' for summation, 'prod' for multiplication or 'kron' for the Kronecker 
# product.
# 
# @return A kernel object
# @export
# @name kernel.comp
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


#' @name composition
#' @rdname composition
#'
#' @title Compositional kernels
#' 
#' @description Construct a new kernel by combining existing kernels 
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

# ~~~~~~~~~~
# basis kernels

#' @title Radial basis function kernel
#' 
#' @description Construct a radial basis function (A.K.A. squared-exponential) kernel.
#' 
#' @details The rbf kernel takes the form:
#' \deqn{k_{rbf} = \sigma^2 exp(-\frac{{(x - x')}^2}{2l^2})}
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
