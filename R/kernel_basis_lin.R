# kernel_lin

#' @title Linear kernel
#' 
#' @description Construct a linear kernel.
#' 
#' @details The linear kernel takes the form:
#' \deqn{k_{lin}(\mathbf{x}, \mathbf{x}') = \sigma_b^2 + \sigma_v^2 (\mathbf{x} - c)(\mathbf{x}' - c)}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active,
#' \eqn{c} determines the value(s) of \eqn{x} through which all realisations pass,
#' \eqn{\sigma_b^2} is a prior over the absolute value of the intercept and
#' \eqn{\sigma_v^2} isa  prior over the slopes of the realisations.
#' 
#' @template kco
#' @template kco_basis
#' @export
#' @name lin
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- lin('temperature')
#' 
#' # and another with two features
#' k2 <- lin(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
lin <- function (columns) {
  
  # construct a linear kernel
  createKernelConstructor('lin',
                          columns,
                          list(sigma2_b = 1,
                               sigma2_v = 1,
                               c = rep(0,
                                       length(columns))),
                          linEval)
  
}

linEval <- function(object, data, newdata = NULL) {
  # evaluate linear kernel against data
  
  # extract from/to data
  data <- getFeatures(object, data, newdata)
  
  x <- data$x
  y <- data$y
  
  # get kernel parameters
  parameters <- object$parameters
  
  # extract lengthscales and variance
  sigma2_b <- parameters$sigma2_b
  sigma2_v <- parameters$sigma2_v
  c <- parameters$c
  
  # subtract the x axis offsets
  x <- sweep(x, 2, c, '-')
  y <- sweep(y, 2, c, '-')
  
  # get distances
  d <- x %*% t(y)
  
  # complete covariance matrix
  covmat <- sigma2_b + sigma2_v * d
  
  # and return
  return (covmat)
  
}