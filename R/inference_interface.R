# interface for fitting Gaussian process models

# this'll need some roxygen
gp <- function(response,
               kernel,
               data,
               family = gaussian,
               mean_function = NULL,
               inducing_data = NULL,
               inference = c('default', 'exact', 'FITC','Laplace')) {
  
  # catch the call
  call <- match.call()
  
  # if there's no mean function, use zeroes
  if (is.null(mean_function)) {
    mean_function <- zeroes
  }
  
  # get the likelihood
  likelihood <- getLikelihood(family)
  
  # match the inference option
  inference <- match.arg(inference)
  
  # get inference function
  inference <- getInference(inference,
                            likelihood)
  
  # check it's a valid kernel
  checkKernel(kernel)  
  
  # do inference
  posterior <- inference(response,
                         data,
                         kernel,
                         likelihood,
                         mean_function,
                         inducing_data)
  
  ans <- list(call = call,
              likelihood = likelihood,
              kernel = kernel,
              mean_function = mean_function,
              data = list(response = response,
                          training_data = data,
                          inducing_data = inducing_data),
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
  
  # if inference is default, look up the default for that likelihood
  if (inference == 'default') {
    inference <- defaultInference(likelihood)
  }
  
  # check inference method is suitable for the likelihood
  checkInference(inference, likelihood)
  
  # fetch the function
  inference <- get(inference)
  
  # return this function
  return (inference)
  
}
