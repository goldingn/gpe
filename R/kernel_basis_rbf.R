# kernel_rbf

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
  createKernelConstructor('rbf',
                          columns,
                          list(sigma = 1,
                               l = rep(1,
                                       length(columns))),
                          rbfEval)
  
}


rbfEval <- function(object, data, newdata = NULL, diag = FALSE) {
  # evaluate rbf kernel against data
  
  # diagonal case
  if (diag) {
    
    # make sure it's symmetric (newdata is null)
    checkSymmetric(newdata)
    
    # if it's fine return sigma squared on the diagonals
    covmat <- diagSigma(object, data)
    
    return (covmat)
    
  }
  
  
  # extract from/to data
  data <- getFeatures(object, data, newdata)
  
  x <- data$x
  y <- data$y
  
  # get kernel parameters
  parameters <- object$parameters
  
  # extract lengthscales and variance
  l <- parameters$l
  sigma <- parameters$sigma
  
  # apply the lengthscale parameters
  x <- sweep(x, 2, l^2, '/')
  y <- sweep(y, 2, l^2, '/')
  
  # get distances
  d <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * exp(-(d ^ 2) / 2)
  
  # and return
  return (covmat)
  
}