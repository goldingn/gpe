# kernel_iid

#' @title IID random effects kernel
#' 
#' @description Construct an iid normal `random effects' kernel.
#' 
#' @param column optionally, a single string giving the name of a \emph{factor} 
#' feature on which the kernel acts. This will be used to access columns from 
#' a dataframe when the kernel is evaluated . If \code{NULL} then the iid 
#' effect acts on each datapoint, rather than on groups.
#'  
#' @details The iid kernel takes the form:
#' \deqn{k_{iid}(\mathbf{x}, \mathbf{x}') = \sigma^2 ( \mathbf{x} \mathbf{x}' ) }
#' where \eqn{\mathbf{x}} are indicator variables with each column containing 
#' 1s for records in a given group and 0s otherwise and \eqn{\sigma^2} is the 
#' overall variance. This is equivalent to a hierarchical or random-effects 
#' model: \eqn{z_j \sim N(0, \sigma^2)} with \eqn{j} indexing the groups of the 
#' specified factor.
#' 
#' In practice, the active column should actually be a single factor and the 
#' indicator variables will be built internally.
#' 
#' @template par_sigma
#' @template kco
#' @export
#' @name iid
#' 
#' @examples
#'  
#' # add a discrete variable to the pressure dataset
#' pressure$group <- factor(sample(letters[1:3], 19, replace = TRUE))
#' 
#' # construct an iid kernel over this
#' k1 <- iid('group')
#' 
#' # evaluate and visualise it
#' image(k1(pressure))
#' 
iid <- function (column = NULL, sigma = 1) {
  
  # throw an error if more than one column is specified
  if (length(column) > 1) {
    stop (paste0('iid kernels can only be constructed on one column at a time',
                 ', perhaps you should try constructing two and summing them?'))
  }
  
  # construct an iid kernel
  createKernelConstructor('iid',
                          column,
                          list(sigma = pos(sigma)),
                          iidEval)
  
}


iidEval <- function(object, data, newdata = NULL, diag = FALSE) {
  
  # evaluate iid kernel against data
  
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
  
  # extract lengthscales and variance
  sigma <- parameters$sigma()
  
  # if columns is null, it's iid on all observations
  if (is.null(object$columns)) {
    
    # get training n
    n <- nrow(data) 
    
    if (is.null(newdata)) {
      
      # if newdata is NULL, it must be the self matrix, so identity
      covmat <- sigma ^ 2 * diag(n)
      
    } else {
      
      # otherwise it's 0s (new data independent so no covariance)
      
      # dimension of evaluation data
      m <- nrow(newdata)
      
      # 0s matrix
      covmat <- matrix(0,
                       nrow = n,
                       ncol = m)
      
    }
    
  } else {
    
    # otherwise, iid on some grouping factor
    
    # extract from/to data, don't vconvert them to matrices
    data <- getFeatures(object, data, newdata, to_matrix = FALSE)
    
    # turn factors into full set of indicator variables
    x <- expandFactor(data$x)
    y <- expandFactor(data$y)
    
    # complete covariance matrix
    covmat <- sigma ^ 2 * x %*% t(y)
    
  }
  
  # return covariance matrix
  return (covmat)
  
}
