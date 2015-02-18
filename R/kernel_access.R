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
  checkKernel(kernel)
  
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
  checkKernel(kernel)
  
  # flatten into a list of the basis kernels
  kernel_list <- flattenKernel(kernel)$data_list
  
  # loop through listing only unique columns
  ans <- unique(unlist(lapply(kernel_list,
                              function(x) x$columns)))
  
  # remove any which are blank
  ans <- ans[ans != ""]
  
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
  checkKernel(kernel)
  
  # extract the type
  ans <- environment(kernel)$object$parameters
  
  # and return it
  return (ans)
  
}

#' @rdname access
#' @param which which of the two \emph{immediate} subkernels of the kernel to extract
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
  checkKernel(kernel)
  type <- getType(kernel) 
  
  # check it's a basis kernel
  checkCompositional(kernel)
  
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
#' @param \dots parameters to be updated and the new values for them to take (see examples).
#' 
#' @export
#' @examples
#' # evaluate and visualise it
#' image(k1(pressure))
#' 
#' # change the length scale
#' k2 <- setParameters(k1, l = 10)
#' getParameters(k2)
#' 
#' # change the lengthhscale and variance
#' k2 <- setParameters(k1, l = 9, sigma = 1.3)
#' getParameters(k2)
#' 
#' # evaluate and visualise the new kernel
#' image(k2(pressure))
#'  
setParameters <- function (kernel, ...) {
  # function to set the values of some or all of the parameters
  # parameter_list must be a named list giving the parameters to be updated
  
  # capture dots argument
  parameter_list <- list(...)
  
  # check the new parameters
  checkParameters(kernel, parameter_list)
  
  # copy over the kernel
  new_kernel <- kernel
  
  # clone the previous environment
  environment(new_kernel) <- environment(kernel)

  # get kernel's object
  object <- getObject(new_kernel)

  # get the existing parameters
  parameters_current <- object$parameters
  
  # loop through each element in parameter_list
  for (i in 1:length(parameter_list)) {
    
    # get the new parameter name
    parameter_name <- names(parameter_list)[i]
    
    # and length
    parameter_len <- length(parameter_list[[i]])
    
    # find the matching parameter
    j <- match(parameter_name,
               names(parameters_current))
    
    # update the parameter
    parameters_current[[j]] <- parameter_list[[i]]
    
  }
  
  # insert parameters back into the kernel
  environment(new_kernel)$object$parameters <- parameters_current
    
  
  # return the kernel
  return (new_kernel)
  
}


#' @rdname access
#' 
#' @param data a dataset from which to calculate the range of values over which to draw example GPs.
#' If \code{NULL}, the range -5 to 5 is used, arbitrarily.
#' @param ndraw the number of (zero-mean) GP realisations from the the kernel to plot.
#' 
#' @export
#' @examples
#'  
#' #plot example GPs from this kernel, applied to the pressure dataset
#' demoKernel(k1, data = pressure)
#' 
demoKernel <- function (kernel, data = NULL, ndraw = 5) {
  
  # plotting example draws from a kernel
  
#   # for now, throw an error the kernel is compositional
#   if (getType(kernel) %in% c('sum', 'prod', 'kron')) {
#     stop ("sorry, this functionality is currently only available for basis kernels")
#   }
  
  # get the columns
  columns <- getColumns(kernel)
  
  D <- length(columns)
  
  if (D != 1) {
    stop ("sorry, only one-dimensional kernels are currently supported")
  }
  
  if (length(columns) == 0) {
    stop ("sorry, only kernels which use a covariate are currently supported")
  }
  
  # need to make this handle discrete data
  
  # get the range of data to plot
  if (is.null(data)) {
    range <- c(-5, 5)
  } else {
    range <- range(data[, columns])
  }
  
  # get some random draws from the kernel
  df_covs <- data.frame(seq(range[1], range[2], len = 1000))
  names(df_covs) <- columns
  
  draws <- rgp(n = ndraw, kernel = kernel, data = df_covs)
  
  plot(draws[, 1] ~ df_covs[[1]],
       type = 'n',
       ylim = range(draws),
       xlab = columns,
       ylab = paste0('f(', columns, ')'))
  
  for (i in 1:ndraw) {
    lines(draws[, i] ~ df_covs[[1]],
          type = 'l',
          lwd = 2,
          col = grey(0.5))
  }
  
  title(main = paste0(ndraw,
                      ' GP realisations from the kernel\n',
                      trimParentheses(parseKernelStructure(kernel))))
  
}
