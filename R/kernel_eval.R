# functions to evaluate kernels

kern.rbf.eval <- function(object, data, newdata = NULL) {
  # evaluate rbf kernel against data
  
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
  covmat <- sigma * exp(-(d ^ 2) / 2)
  
  # and return
  return (covmat)
  
}

kern.lin.eval <- function(object, data, newdata = NULL) {
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

kern.per.eval <- function(object, data, newdata = NULL) {
  # evaluate rbf kernel against data
  
  # extract from/to data
  data <- getFeatures(object, data, newdata)
  
  x <- data$x
  y <- data$y
  
  # get kernel parameters
  parameters <- object$parameters
  
  # extract lengthscales and variance
  p <- parameters$p
  l <- parameters$l
  sigma2 <- parameters$sigma2
  
  # get distances
  d <- fields::rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma2 * exp(-(2 * sin(pi * d / p) ^ 2) / l ^ 2)
  
  # and return
  return (covmat)
  
}

kern.iid.eval <- function(object, data, newdata = NULL) {
  
  # evaluate iid kernel against data
    
  # get kernel parameters
  parameters <- object$parameters
  
  # extract lengthscales and variance
  sigma <- parameters$sigma
  
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
  
    # extract from/to data
    data <- getFeatures(object, data, newdata)
    
    # turn factors into full set of indicator variables
    x <- expandFactor(data$x)
    y <- expandFactor(data$y)
        
    # complete covariance matrix
    covmat <- sigma ^ 2 * x %*% t(y)
  
  }
  
  # return covariance matrix
  return (covmat)
  
}

kern.comp.eval <- function(object, data, newdata, operation) {
  
  stopifnot(operation %in% c('sum', 'prod', 'kron'))
  
  # extract the sub kernels
  kernel1 <- getObject(object$kernel1)
  kernel2 <- getObject(object$kernel2)
  
  # evaluate them
  covmat1 <- evalKernel(kernel1, data, newdata)
  covmat2 <- evalKernel(kernel2, data, newdata)
  
  # get the compositional covariance matrix
  if (operation == 'sum') {
    
    covmat <- covmat1 + covmat2
    
  } else if (operation == 'prod') {
    
    covmat <- covmat1 * covmat2
    
  } else {
    
    covmat <- covmat1 %x% covmat2
    
  }
  
  # return the covariance matrix
  return (covmat)
}


# evaluation master function
evalKernel <- function (object, data, newdata = NULL) {
  
  # evaluate the kernel against data
  
  
  # execute the different cases
  
  covmat <- switch(object$type,
                   sum = kern.comp.eval(object,
                                        data,
                                        newdata,
                                        'sum'),
                   prod = kern.comp.eval(object,
                                         data,
                                         newdata,
                                         'prod'),
                   kron = kern.comp.eval(object,
                                         data,
                                         newdata,
                                         'kron'),
                   rbf = kern.rbf.eval(object,
                                       data,
                                       newdata),
                   lin = kern.lin.eval(object,
                                       data,
                                       newdata),
                   per = kern.per.eval(object,
                                       data,
                                       newdata),
                   iid = kern.iid.eval(object,
                                       data,
                                       newdata))
}
