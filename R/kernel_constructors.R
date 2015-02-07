# kernel constructor functions

# compositional kernel function - should this be exported?
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

rbf <- function (columns) {
  
  # construct an rbf kernel
  
  # create the model object and initialize parameters
  object <- list(type = 'rbf',
                 columns = columns,
                 parameters = list(sigma = 1,
                                   l = rep(1,
                                           length(columns))))
  
  # create a function to return
  ans <- function(data, newdata = NULL) {
    
    evalKernel(object, data, newdata)    
    
  }
  
  # tell this function it's now a kernel
  ans <- as.kernel(ans)
  
  return(ans)
  
}

