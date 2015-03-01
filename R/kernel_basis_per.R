# kernel_per

#' @title Periodic kernel
#' 
#' @description Construct a periodic kernel.
#' 
#' @details The periodic kernel takes the form:
#' \deqn{k_{per}(\mathbf{x}, \mathbf{x}') = \sigma^2 exp \left(-\frac{2sin^2(\pi | \mathbf{x} - \mathbf{x}' | /p)}{l^2} \right)}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active,
#' \eqn{p} determines the periodicity (distance between successive peaks),
#' \eqn{l} is a characteristic lengthscale, as in the rbf kernel, and \eqn{\sigma^2}
#' is the amplitude of the signal
#' 
#' @template kco
#' @template kco_basis
#' @export
#' @name per
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- per('temperature')
#' 
#' # and another with two features
#' k2 <- per(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
per <- function (columns) {
  
  # construct a periodic kernel
  createKernelConstructor('per',
                          columns,
                          list(p = pos(1),
                               l = pos(1),
                               sigma = pos(1)),
                          perEval)
  
}


perEval <- function(object, data, newdata = NULL, diag = FALSE) {
  # evaluate periodic kernel against data
  
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
  p <- parameters$p()
  l <- parameters$l()
  sigma <- parameters$sigma()
  
  # get distances
  d <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * exp(-(2 * sin(pi * d / p) ^ 2) / l ^ 2)
  
  # and return
  return (covmat)
  
}