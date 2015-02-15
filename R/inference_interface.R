# interface for fitting Gaussian process models

# this'll need some roxygen
gp <- function(response,
               kernel,
               data,
               family = gaussian,
               mean_function = NULL,
               inducing_data = NULL,
               inference = c('default', 'exact', 'FITC')) {

  # capture family as a string
  family <- deparse(substitute(family))
  
  # if there's no mean function, use zeroes
  if (is.null(mean_function)) {
    
    # return the default (zeroes)
    mean_function <- zeroes
    
  }
  
  # get the inference
  inference <- match.arg(inference)

  # if it's default, look up the default
  if (inference == 'default') {
    inference <- defaultInference(family)
  }
    
  # get inference function
  inference_fun <- getInference(inference, family)

  # check it's a valid kernel
  checkKernel(kernel)  
  
  # do inference
  posterior <- inference_fun(y,
                   data,
                   kernel,
                   mean_function,
                   inducing_data)
  
  ans <- list(family = family,
              posterior = posterior)
  
  class(ans) <- 'gp'
  
  return (ans)
  
}

checkInference <- function (inference, family) {
  
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

# return an inference function, given an inference method and a family
getInference <- function (inference, family) {
 
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
