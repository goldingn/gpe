# kernel_rq

#' @title Rational quadratic kernel
#' 
#' @description Construct a rational quadratic kernel.
#' 
#' @details The rational quadratic kernel takes the form:
#' \deqn{k_{rq}(\mathbf{x}, \mathbf{x}') =
#'  \sigma^2 (\left[
#'    1 + {\sum\limits_{d=1}^D \left(
#'        \frac{(x_d - x_d')}{2 \alpha l_d^2}
#'        \right)}^2
#'        \right] ^ -\alpha}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active, 
#' \eqn{l_d} are the characteristic lengthscales for each covariate (column) 
#' \eqn{x_d}, \eqn{\sigma^2} is the overall variance and \eqn{\alpha} controls the
#' roughness.
#' 
#' Larger values of \eqn{l_i} correspond to functions in which change less 
#' rapidly over the values of the covariates.
#' 
#' @template par_sigma
#' @template par_l
#' @template par_alpha
#' @template kco
#' @template kco_basis
#' @export
#' @name rq
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- rq('temperature')
#' 
#' # and another with two features
#' k2 <- rq(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
rq <- function (columns, sigma = 1, l = rep(1, length(columns)), alpha = 1) {
  
  # construct an rbf kernel
  createKernelConstructor('rq',
                          columns,
                          list(sigma = pos(sigma),
                               l = pos(l),
                               alpha = pos(alpha)),
                          rqEval)
  
}


rqEval <- function(object, data, newdata = NULL, diag = FALSE) {
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
  l <- parameters$l()
  sigma <- parameters$sigma()
  alpha <- parameters$alpha()
  
  # apply the lengthscale parameters
  x <- sweep(x, 2, 2 * alpha * l ^ 2, '/')
  y <- sweep(y, 2, 2 * alpha * l ^ 2, '/')
  
  # get distances
  d <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * (1 + (d ^ 2)) ^ -alpha
  
  # and return
  return (covmat)
  
}