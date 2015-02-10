# kernel constructor functions

#' @name composition
#' @rdname composition
#'
#' @title Compositional kernels
#' 
#' @description Construct a new kernel by combining existing kernels, either 
#' by summation, multiplication or the kronecker product.
#' 
#' Summation and multiplication require that the covariance matrices produced
#' the two kernels have the same dimension (same number of rows in the input
#' columns) and result in a kernel with the same dimension as its inputs.
#' 
#' The kronecker product doesn't require that the input functions have the 
#' same dimension, and the dimension of the output is the product of the 
#' dimensions of the inputs (i.e. an m-by-m matrix kroneckered with an 
#' n-by-n matrix gives rise to an nm-by-nm matrix).
#' 
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

# underlying compositional kernel function
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


# ~~~~~~~~~~
# basis kernels

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
  
  # create the model object and initialize parameters
  object <- list(type = 'lin',
                 columns = columns,
                 parameters = list(sigma2_b = 1,
                                   sigma2_v = 1,
                                   c = rep(0,
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

#' @title Periodic kernel
#' 
#' @description Construct a periodic kernel.
#' 
#' @details The periodic kernel takes the form:
#' \deqn{k_{lin}(\mathbf{x}, \mathbf{x}') = \sigma^2 exp \left(-\frac{2sin^2(\pi | \mathbf{x} - \mathbf{x}' | /p)}{l^2} \right)}
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
  
  # create the model object and initialize parameters
  object <- list(type = 'per',
                 columns = columns,
                 parameters = list(p = 1,
                                   l = 1,
                                   sigma2 = 1))
  
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


#' @title IID random effects kernel
#' 
#' @description Construct an iid normal 'random effect' kernel.
#' 
#' @param column optionally, a single string giving the name of a \emph{factor} 
#' feature on which the kernel acts. This will be used to access columns from 
#' a dataframe when the kernel is evaluated . If \code{NULL} then the iid 
#' effect acts on each datapoint, rather than on groups.
#'  
#' @details The iid kernel takes the form:
#' \deqn{k_{iid}(\mathbf{x}, \mathbf{x}') = \sigma^2 \mathbf{x} \mathbf{x}'}
#' where \eqn{\mathbf{x}} are indicator variables with each column containing 
#' 1s for records in a given group and 0s otherwise and \eqn{\sigma^2} is the 
#' overall variance. This is equivalent to a hierarchical or random-effects 
#' model: \eqn{z_j ~ N(0, \sigma2)} with \eqn{j} indexing the groups of the 
#' specified factor.
#' 
#' In practice, the active column should actually be a single factor and the 
#' indicator variables will be built internally.
#' 
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
iid <- function (column = NULL) {
  
  # construct an iid kernel

  # throw an error if more than one column is specified
  if (length(column) > 1) {
    stop (paste0('iid kernels can only be constructed on one column at a time',
          ', perhaps you should try constructing two and summing them?'))
  }
  
  # create the model object and initialize parameters
  object <- list(type = 'iid',
                 columns = column,
                 parameters = list(sigma = 1))
  
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