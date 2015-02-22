# interface for fitting Gaussian process models

# this'll need some roxygen
gp <- function(response,
               kernel,
               data,
               family = gaussian,
               mean_function = NULL,
               inducing_data = NULL,
               inference = c('default', 'exact', 'FITC')) {


  # need to decide on how family is specified - use glm family and
  # just pull out likelihood & link?
  # need a binomial wrapper around bernoulli likelihood
  
  # if there's no mean function, use zeroes
  if (is.null(mean_function)) {
    
    # return the default (zeroes)
    mean_function <- zeroes
    
  }
  
  # get the likelihood
  likelihood <- getLikelihood(family)
  
  # get the inference
  inference <- match.arg(inference)

  # if it's default, look up the default
  if (inference == 'default') {
    inference <- defaultInference(likelihood)
  }
    
  # get inference function
  inference_fun <- getInference(inference, likelihood)

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

checkInference <- function (inference, likelihood) {
  
  if (inference %in% c('exact', 'FITC') &
        likelihood$name != 'likelihood_gaussian_identity') {
    stop ('exact and FITC inference are only possible with the gaussian family and identity link')
  }
  
}

defaultInference <- function (likelihood) {
  
  inference <- switch(likelihood$name,
                      likelihood_gaussian_identity = 'inference_direct_exact',
                      likelihood_binomial_logit = 'inference_laplace_exact',
                      likelihood_binomial_probit = 'inference_laplace_exact',)
  
  # if nothing found, throw an error
  if (is.null(inference)) {
    stop (paste0('no default inference method found for likelihood ',
                 likelihood$name))
  }
  
  return (inference)
  
}

# return an inference function, given an inference method and a family
getInference <- function (inference, likelihood) {
 
  # check inference method is suitable for the likelihood
  checkInference(inference, likelihood)
  
  # fetch the function
  inference <- get(inference)
    
  # return this function
  return (inference)

}

