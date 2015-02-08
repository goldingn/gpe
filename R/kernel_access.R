# Funtions to access parts of kernels

#' @name access
#' @rdname access
#'
#' @title Interaction with kernel objects
#' 
#' @description Access and interact with kernel objects. 
#' These functions allow you to report the type of kernel (\code{getType}), 
#' view the overall kernel structure (\code{?}), extract subkernels (\code{getSubKernel})  
#' report kernel parameters (\code{getParameters}) and update them (\code{setParameters}).
#' 
#' @template kac_kernel
#' 
#' @examples
#'  
#' # construct a kernel with one feature
#' k1 <- rbf('temperature')
#'  
NULL

#' @rdname access
#' @export
#' @examples
#'  
#' # get the kernel type
#' getType(k1)
#'  
getType <- function (kernel) {
  
  # get the type of kernel
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$type
  
  # and return it
  return (ans)
  
}
#' @rdname access
#' @export
#' @examples
#'  
#' # get the names of the columns the kernel acts on
#' getColumns(k1)
#'  
getColumns <- function (kernel) {
  
  # get the columns used to construct the kernel
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$columns
  
  # and return it
  return (ans)
  
}

#' @rdname access
#' @export
#' @examples
#'  
#' # get the parameters of the kernel
#' params <- getParameters(k1)
#' params
#'  
getParameters <- function (kernel) {
  
  # get the kernel parameters
  
  # check it's a kernel
  stopifnot(is.kernel(kernel))
  
  # extract the type
  ans <- environment(kernel)$object$parameters
  
  # and return it
  return (ans)
  
}

#' @rdname access
#' @param which Which of the two \emph{immediate} subkernels of the kernel to extract
#' @export
#' @examples
#'  
#' # build a compositional kernel
#' k2 <- k1 + k1 * k1
#'  
#' # extract a subkernel
#' k3 <- getSubKernel(k2, 1)
#'  
getSubKernel <- function (kernel, which = 1) {
  
  # extract the subkernels of compositional kernels
  
  # kernel must be a sum, prod or kron kernel and which can be either 1 or 2
  
  # check the class and type
  stopifnot(is.kernel(kernel))
  type <- getType(kernel) 
  
  # if it's not a compositional kernel give a nice error
  if (!(type %in% c('sum', 'prod', 'kron'))) {
    
    stop(paste0('only compositional kernel have sub-kernels,',
                'this is a basis kernel of type: ',
                type))
    
  }
  
  # check the subkernel selection
  if (!(which %in% 1:2)) {
    stop (paste0('this function can only access the *immediate* sub-kernels',
                 ' of a compositional kernel, so which can only be 1 or 2.',
                 ' You entered ',
                 which))
  }
  
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

#' @rdname access
#' 
#' @param parameter_list A named list with names giving the parameters to be 
#' updated and values giving the new values they should take on.
#' 
#' @export
#' @examples
#' # evaluate and visualise it
#' image(k1(pressure))
#' 
#' # change the length scale
#' params$l <- 100
#' k2 <- setParameters(k1, params)
#' 
#' # evaluate and visualise the new kernel
#' image(k1(pressure))
#'  
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
      if (!(parameter_name %in% names(parameters_current))) {
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
