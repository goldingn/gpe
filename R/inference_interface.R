# interface for fitting Gaussian process models

#' @title Gaussian process models
#' 
#' @description Fiting (latent) Gaussian process models using a range of 
#' inference methods.
#' 
#' @param response A vector or matrix containing values of the response 
#' variable to be modelled.
#'
#' @param kernel A kernel object.
#' 
#' @param data A data frame containing the covariates against which to model 
#' the response variable. This must have the same number of rows as 
#' \code{response} and contain named variables matching those referred to by
#' \code{kernel}.
#' 
#' @param family A \code{\link{family}} object giving the likelihood and link 
#' function to be used to fit the model. Currently only \code{gaussian} is 
#' supported.
#' 
#' @param mean_function An optional function specifying the prior over the mean
#' of the gp, in other words a 'first guess' at what the true function is.
#' This must act on a dataframe with named variables matching some of those in 
#' \code{data} and return a vector giving a single value for each row in the 
#' dataframe. Note that this function must return a prediction on the scale of 
#' the link, rather than the response. If \code{NULL} then a prior mean is 
#' assumed to be 0 for all observations.
#' 
#' @param inducing_data An optional dataframe containing the locations of 
#' inducing points to be used when carrying out sparse inference (e.g. FITC).
#' This must contain variables with names matching those referenced by 
#' \code{kernel} and \code{mean_function}. This should have fewer rows than 
#' \code{data} and \code{response}. 
#' 
#' @param inference A string specifying the inference method to be used to 
#' estimate the values of the latent parameters. If \code{'default'} an 
#' appropriate method is picked for the likelihood specified. See details 
#' section for the list of default inference methods.
#' 
#' @return A fitted gp object for which there aren't yet any associated 
#' functions. But there will be.
#' 
#' @details The default inference method for a model with the family 
#' \code{gaussian(link = 'identity')} is full, direct inference.
#' 
#' @export
#' @name gp
#' 
#' @examples
#' 
#' # make some fake data
#' n <- 100  # observations
#' m <- 10  # inducing points
#' 
#' # dataframes
#' df <- data.frame(x = sort(runif(n, -5, 5)))
#' inducing_df <- data.frame(x = sort(runif(m, -5, 5)))
#' prediction_df <- data.frame(x = seq(min(df$x), max(df$x), len = 500))
#' 
#' # fake Gaussian response data
#' f <- sin(df$x)
#' y <- rnorm(n, f, 1)
#' 
#' # construct a kernel (including the observation error)
#' kernel <- rbf('x') + iid()
#' 
#' # fit a full (non-sparse) GP model (without updating the hyperparameters) 
#' # as this is the default
#' m1 <- gp(y, kernel, df, gaussian)
#' 
#' # fit another with FITC sparsity
#' m2 <- gp(y, kernel, df, gaussian, inference = 'FITC', inducing_data = df)#inducing_df)
#' # summary stats, other associated functions still to come
#' 
# this'll need some roxygen
gp <- function(response,
               kernel,
               data,
               family = gaussian,
               mean_function = NULL,
               inducing_data = NULL,
               inference = c('default', 'full', 'FITC','Laplace'),
               verbose = FALSE) {
  
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
                            likelihood,
                            inducing_data)
  
  # check it's a valid kernel
  checkKernel(kernel)  
  
  # do inference
  posterior <- inference(response,
                         data,
                         kernel,
                         likelihood,
                         mean_function,
                         inducing_data,
                         verbose = verbose)
  
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

checkInference <- function (inference, likelihood, inducing_data) {
  
  if (inference %in% c('full', 'FITC') &
        likelihood$name != 'likelihood_gaussian_identity') {
    stop ('full and FITC inference are only possible with the gaussian family and identity link')
  }
  
  if (inference %in% c('FITC', 'LaplaceFITC') &
        is.null(inducing_data)) {
    stop ('If FITC inference is required, inducing_data must be provided.')
  }
  
}

defaultInference <- function (likelihood) {
  
  inference <- switch(likelihood$name,
                      likelihood_gaussian_identity = 'full',
                      likelihood_binomial_logit = 'Laplace',
                      likelihood_binomial_probit = 'Laplace')
  
  # if nothing found, throw an error
  if (is.null(inference)) {
    stop (paste0('no default inference method found for likelihood ',
                 likelihood$name))
  }
  
  return (inference)
  
}

# return an inference function, given an inference method and a family
getInference <- function (inference, likelihood, inducing_data) {
  
  # if inference is default, look up the default for that likelihood
  if (inference == 'default') {
    inference <- defaultInference(likelihood)
  }
  
  # check inference method is suitable for the likelihood
  checkInference(inference, likelihood, inducing_data)
  
  # get the inference name
  inference_name <- switch(inference,
                           FITC = 'direct_fitc',
                           full = 'direct_full',
                           LaplaceFITC = 'laplace_fitc',
                           Laplace = 'laplace_full')
  
  # if the inference method wasn't found throw an error
  if (is.null(inference_name)) {
    stop (paste0('Inference method ',
                 inference,
                 ' not found.'))
  }
  
  # prepend 'inference' to the name
  inference_name <- paste0('inference_',
                          inference_name)
  
  # fetch the function
  inference <- get(inference_name)
  
  # return this function
  return (inference)
  
}
