gp <- function(response,
               kernel,
               data,
               family = gaussian,
               mean_function = NULL,
               inducing_points = NULL,
               inference = c('default', 'exact', 'FITC')) {

  # capture family as a string
  family <- deparse(substitute(family))
  
  # and inference
  inference <- match.arg(inference)
  # check it's a valid kernel
  checkKernel(kernel)  
  
  # get inference function
  inference <- getInference(inference, family)

  # do inference
  ans <- inference(y,
                   data,
                   inducing_data,
                   kernel,
                   mean_function)
  
  class(ans) <- 'gp'
  
  return (ans)
  
}

checkInference(inference, family) {
  
  if (inference %in% c('exact', 'FITC') &
        family != 'gaussian') {
    stop ('exact and FITC inference are only possible with the gaussian family')
  }
  
}

defaultInference <- function (family) {
  
  inference <- switch(family,
                      gaussian = 'exact')
  
  # if nothing found, throw an error
  if (is.null(inference)) {
    stop (paste0('no default inference method found for family: ',
                 family))
  }
  
  return (inference)
  
}

getInference <- function (inference, family) {

  # if it's default, look up the default
  if (inference == 'default') {
    inference <- defaultInference(family)
  }
  
  # check inference method is suitable for the family
  checkInference(inference, family)
  
  # switch inference to method
  inference <- switch(inference,
                      exact = infExact,
                      FITC = infFITC)
  
  if (is.null(inference)) {
    stop ('inference method not available')
  }
  
  return (inference)

}