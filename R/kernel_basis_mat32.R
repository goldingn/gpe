# kernel_mat32

#' @title Matern 3/2 kernel
#' 
#' @description Construct a Matern 3/2 kernel.
#' 
#' @details The Matern 3/2 kernel takes the form:
#' \deqn{k_{mat32}(\mathbf{x}, \mathbf{x}') = \sigma^2 (1 + \sqrt{3} \mathbf{r}) exp(-\sqrt{3} \mathbf{r})} 
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
#' @name mat32
#' 
#' @examples
#' # construct a kernel with one feature
#' k1 <- mat32('temperature')
#' 
#' # and another with two features
#' k2 <- mat32(c('temperature', 'pressure'))
#' 
#' # evaluate them on the pressure dataset
#' image(k1(pressure))
#' image(k2(pressure))
#' 
mat32 <- function (columns, sigma = 1, l = rep(1, length(columns))) {
  
  # construct an rbf kernel
  createKernelConstructor('mat32',
                          columns,
                          list(sigma = pos(sigma),
                               l = pos(l)),
                          mat32Eval)
  
}


mat32Eval <- function(object, data, newdata = NULL, diag = FALSE) {
  # evaluate Matern 3/2 kernel against data
  
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
  x <- sweep(x, 2, l ^ 2, '/')
  y <- sweep(y, 2, l ^ 2, '/')
  
  # get distances
  r <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma ^ 2 * (1 + sqrt(3) * r) * exp(-sqrt(3) * r)
  
  # and return
  return (covmat)
  
}
