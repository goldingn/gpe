# functions to evaluate and access parts of kernels

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
  y <- sweep(y, 2, l, '/')
  
  # get distances
  d <- rdist(x, y)
  
  # complete covariance matrix
  covmat <- sigma * exp(-(d ^ 2) / 2)
  
  # and return
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
                                       newdata))
  
  return (covmat)
  
}


# functions to access bits of kernels

getType <- function (kernel) {
  
  # get the type of kernel
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$type
  
  # and return it
  return (ans)
  
}

getColumns <- function (kernel) {
  
  # get the columns used to construct the kernel
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$columns
  
  # and return it
  return (ans)
  
}

getParameters <- function (kernel) {
  
  # get the kernel parameters
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$parameters
  
  # and return it
  return (ans)
  
}

getSubKernel <- function (kernel, which = 1) {
  
  # extract the subkernels of compositional kernels
  
  # kernel must be a sum, prod or kron kernel and which can be either 1 or 2
  
  # check the class and type
  stopifnot(is.kernel(kernel))
  stopifnot(getType(kernel) %in% c('sum', 'prod', 'kron'))
  
  # check the subkernel selection
  stopifnot(which %in% 1:2)
  
  # extract the chosen subkernel
  if (which == 1) {
    kernel <- environment(kernel)$object$kernel1
  } else {
    kernel <- environment(kernel)$object$kernel2
  }
  
  # return the kernel
  return (kernel)
  
}

getObject <- function (kernel) {
  
  # get a kernel object, from a kernel
  ans <- environment(kernel)$object
  
  # return this
  return (ans)
  
}

getFeatures <- function (object, data, newdata) {
  
  # extract and tidy data for kernel evaluation
  
  # if no evaluation data provided, use the training data
  if (is.null(newdata)) newdata <- data
  
  # get the kernel's active columns
  columns <- object$columns
  
  x <- data[, columns]
  y <- newdata[, columns]
  
  return(list(x = as.matrix(x),
              y = as.matrix(y)))
  
}

setParameters <- function (kernel, parameter_list) {
  # function to set the values of some or all of the parameters
  # parameter_list must be a named list giving the parameters to be updated
  
  # copy over the kernel
  new_kernel <- kernel
  
  # clone the previous environment
  # This will clone e1
  environment(new_kernel) <- as.environment(as.list(environment(kernel),
                                                    all.names = TRUE))
  
  # get kernel's object
  object <- getObject(new_kernel)
  
  # first check the kernel has parameters to set
  if (object$type %in% c('prod', 'sum', 'kron')) {
    
    stop ('compositional kernels have no parameters to set')
    
  } else {
    
    # otherwise, get the existing parameters
    parameters_current <- object$parameters
    
    # check each element in parameter_list
    for (i in 1: length(parameter_list)) {
      
      # get the new parameter name
      parameter_name <- names(parameter_list)[i]
      
      # and length
      parameter_len <- length(parameter_list[[i]])
      
      # check it matches a valid kernel parameter 
      if (!(names(parameter_list) %in% names(parameters_current))) {
        stop (paste0(parameter_name,
                     'is not a valid parameter.\nValide parameters are:',
                     names(parameters_current)))
      }
      
      # otherwise find matching parameter
      j <- match(parameter_name,
                 names(parameters_current))
      
      # get length of this target parameter
      target_parameter_len <- length(parameters_current[[j]])
      
      # check it is of the correct dimension
      if (target_parameter_len != parameter_len) {
        stop (paste0(parameter_name,
                     'is of length ',
                     parameter_len,
                     'but should be of length ',
                     target_parameter_len))
      }
      
      # --------------------------------
      # check support of the parameters!
      # --------------------------------
      
      # otherwise, update the parameters
      parameters_current[[j]] <- parameter_list[[i]]
      
    }
    
    # insert parameters back into the kernel
    environment(new_kernel)$object$parameters <- parameters_current
    
  }
  
  # return the kernel
  return (new_kernel)
  
}
