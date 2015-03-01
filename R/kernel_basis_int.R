# kernel_basis_int

#' @title Intercept kernel
#' 
#' @description Construct an intercept kernel
#' 
#' @details The intercept kernel takes the form:
#' \deqn{k_{int}(\mathbf{x}, \mathbf{x}') = \sigma^2}
#' where \eqn{\sigma^2} is the variance of the prior over the value of the 
#' intercept. This is equivalent to a normal prior over the value of an 
#' intercept term: \eqn{\alpha ~ N(0, \sigma^2)}.
#' 
#' @template kco
#' @export
#' @name int
#' 
#' @examples
#'  
#' # construct an intercept kernel
#' k1 <- int()
#' 
#' # evaluate and visualise it
#' image(k1(pressure))
#' 
int <- function () {
  
  # construct an iid kernel
  createKernelConstructor('int',
                          '',
                          list(sigma = pos(1)),
                          intEval)
  
}


intEval <- function(object, data, newdata = NULL, diag = FALSE) {
  
  # evaluate intercept kernel against data
  
  # diagonal case
  if (diag) {
    
    # make sure it's symmetric (newdata is null)
    checkSymmetric(newdata)
    
    # if it's fine return sigma squared on the diagonals
    covmat <- diagSigma(object, data)
    
    return (covmat)
    
  }
  # get kernel parameters
  parameters <- object$parameters
  
  # extract variance
  sigma <- parameters$sigma()
  
  # get training n
  n_x <- nrow(data) 

  # and prediction n
  if (is.null(newdata)) {
    n_y <- n_x
  } else {
    n_y <- nrow(newdata)
  }
      
  covmat <- sigma ^ 2 * matrix(1,
                       nrow = n_x,
                       ncol = n_y)
      
  # return covariance matrix
  return (covmat)
  
}
