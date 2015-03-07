# Funtions to access parts of kernels

#' @name access
#' @rdname access
#'
#' @title Interaction with kernel objects
#' 
#' @description Access and interact with kernel objects. 
#' These functions allow you to report the type of kernel (\code{getType}), 
#' extract subkernels (\code{getSubKernel}) report kernel parameters 
#' (\code{getParameters}) and update them (\code{setParameters}).
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
#' 
#' @param continuous whether the parameter values should be reported (and set)
#' on the scale of their continuous transformation rather than their true value
#' (the default).

#' @export
#' @examples
#'  
#' # get the parameters of the kernel
#' params <- getParameters(k1)
#' params
#'  
getParameters <- function (kernel, continuous = FALSE) {
  
  # get the kernel parameters
  
  # check it's a kernel
  checkKernel(kernel)
  
  # extract the parameter objects
  ans <- environment(kernel)$object$parameters
  
  # evaluate them all
  ans <- lapply(ans,
                function(x) x(continuous = continuous))
  
  # and return them
  return (ans)
  
}

# 
getObject <- function (kernel) {
  
  # check it first
  checkKernel(kernel)
  
  # get a kernel object, from a kernel
  ans <- environment(kernel)$object
  
  # return this
  return (ans)
  
}

# given a kernel object (such as one produced by getObject)
# re-create the kernel from which it came
# The kernel can be either compositional or basis.
setObject <- function (object) {
  
  is_comp <- all(names(object) == c('type',
                                    'kernel1',
                                    'kernel2'))
  
  is_basis <- all(names(object) == c('type',
                                     'columns',
                                     'parameters'))
  
  if ((!is_comp & !is_basis) | (is_comp & is_basis)) {
    stop ('apparently not a valid kernel object')
  }
  
  # if it's compositional, just recreate it
  if (is_comp) {
    kernel <- get(object$type)(object$kernel1, object$kernel2)
  } else {
    # otherwise, get the type and columns
    kernel <- get(object$type)(object$columns)
    
    # get the parameter values
    params <- lapply(object$parameters,
                     function(x) x())
    
    # update the parameters with these
    kernel <- do.call(setParameters, c(kernel, params))
  }
  
  # check it's valid
  checkKernel(kernel)
  
  return (kernel)
  
}

getFeatures <- function (object, data, newdata, to_matrix = TRUE) {
  
  # extract and tidy data for kernel evaluation,
  # optionally (default) convert response to a matrix
  
  # if no evaluation data provided, use the training data
  if (is.null(newdata)) newdata <- data
  
  # get the kernel's active columns
  columns <- object$columns
  
  x <- data[, columns]
  y <- newdata[, columns]
  
  # optionally convert to matrices
  if (to_matrix) {
    x <- as.matrix(x)
    y <- as.matrix(y)
  }
  
  return(list(x = x,
              y = y))
  
} 

#' @rdname access
#' 
#' @param \dots parameters to be updated and the new values for them to take 
#' (see examples). Note that these must be passed as values, not as parameter
#' objects
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
setParameters <- function (kernel, ..., continuous = FALSE) {
  # function to set the values of some or all of the parameters
  
  # capture dots argument
  parameter_list <- list(...)
  
  # check the new parameters
  checkParameters(kernel, parameter_list)
  
  # get kernel's object
  object <- getObject(kernel)
  
  # get the existing parameters
  parameters_current <- object$parameters
  
  # convert these to their values
  parameters_current_values <- lapply(parameters_current,
                                      function(x, continuous) x(continuous),
                                      continuous)
  
  # loop through each element in parameter_list
  for (i in 1:length(parameter_list)) {
    
    # get the new parameter name
    parameter_name <- names(parameter_list)[i]
    
    # and length
    parameter_len <- length(parameter_list[[i]])
    
    # find the matching parameter
    j <- match(parameter_name,
               names(parameters_current_values))
    
    # update the parameter
    parameters_current_values[[j]] <- parameter_list[[i]]
    
  }
  
  # create a new kernel with these parameter values
  new_kernel <- do.call(getType(kernel),
                        c(list(columns = getColumns(kernel)),
                               parameters_current_values))
  
  # return the kernel
  return (new_kernel)
  
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
  
  # define colour palette
  col <- RColorBrewer::brewer.pal(9, 'Set1')
  
  # repeat it multiple times if necessary
  if (ndraw > 9) {
    col <- rep(col, ceiling(ndraw / 9))
  }
  
  plot(draws[, 1] ~ df_covs[[1]],
       type = 'n',
       ylim = range(draws),
       xlab = columns,
       ylab = paste0('f(', columns, ')'))
  
  for (i in 1:ndraw) {
    lines(draws[, i] ~ df_covs[[1]],
          type = 'l',
          lwd = 2,
          col = col[i])
  }
  
  title(main = paste0(ndraw,
                      ' GP realisations from the kernel\n',
                      trimParentheses(parseKernelStructure(kernel))))
  
}
