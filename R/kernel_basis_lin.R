# kernel_lin

#' @title Linear kernel
#' 
#' @description Construct a linear kernel.
#' 
#' @details The linear kernel takes the form:
#' \deqn{k_{lin}(\mathbf{x}, \mathbf{x}') = \sigma^2 (\mathbf{x} - c)(\mathbf{x}' - c)}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active,
#' \eqn{c} determines the value(s) of \eqn{x} through which all realisations 
#' pass and \eqn{\sigma^2} is a prior over the slopes of the realisations.
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
                          list(sigma = 1,
                               c = rep(0,
                                       length(columns))),
                          linEval)
  
}

linEval <- function(object, data, newdata = NULL, diag = FALSE) {
  # evaluate linear kernel against data
  
  # diagonal case
  if (diag) {
    
    # make sure it's symmetric (newdata is null)
    checkSymmetric(newdata)
    
    # throw an error as it isn't implemented yet
    stop ('diag not implemented for linear kernel yet')
    
  }
  
  # extract from/to data
  data <- getFeatures(object, data, newdata)
  
  x <- data$x
  y <- data$y
  
  # get kernel parameters
  parameters <- object$parameters
  
  # extract lengthscales and variance
  sigma <- parameters$sigma
  c <- parameters$c
  
  # subtract the x axis offsets
  x <- sweep(x, 2, c, '-')
  y <- sweep(y, 2, c, '-')
  
  # get distances
  d <- x %*% t(y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * d
  
  # and return
  return (covmat)
  
}