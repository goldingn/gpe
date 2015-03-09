# kernel_mat52

#' @title Matern 5/2 kernel
#' 
#' @description Construct a Matern 5/2 kernel.
#' 
#' @details The Matern 5/2 kernel takes the form:
#' \deqn{k_{mat52}(\mathbf{x}, \mathbf{x}') = \sigma^2  \left(1 + \sqrt{5} \mathbf{r} + \frac{5 \mathbf{r} ^ 2}{3}\right) exp(- \sqrt{5} \mathbf{r})} 
#' \deqn{\mathbf{r} = {\sqrt{\sum\limits_{d=1}^D \left(\frac{(x_d - x_d')}{2l_d^2}\right) ^ 2}}}
#' where \eqn{\mathbf{x}} are the covariates on which the kernel is active, \eqn{l_d} 
#' are the characteristic lengthscales for each covariate (column) \eqn{x_d} 
#' and \eqn{\sigma^2} is the overall variance.
#' Larger values of \eqn{l_i} correspond to functions which change less 
#' rapidly over the values of the covariates.
#' 
#' @template par_sigma
#' @template par_l
#' @template kco
#' @template kco_basis
#' @export
#' @name mat52
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- mat52('temperature')
#' 
#' # and another with two features
#' k2 <- mat52(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
mat52 <- function (columns, sigma = 1, l = rep(1, length(columns))) {
  
  # construct an rbf kernel
  createKernelConstructor('mat52',
                          columns,
                          list(sigma = pos(sigma),
                               l = pos(l)),
                          mat52Eval)
  
}


mat52Eval <- function(object, data, newdata = NULL, diag = FALSE) {
  # evaluate Matern 5/2 kernel against data
  
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
  
  # apply the lengthscale parameters
  x <- sweep(x, 2, 2 * l ^ 2, '/')
  y <- sweep(y, 2, 2 * l ^ 2, '/')
  
  # get distances
  r <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * (1 + sqrt(5) * r + 5 * (r ^ 2) / 3) * exp(-sqrt(5) * r)
  
  # and return
  return (covmat)
  
}
