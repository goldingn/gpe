# gp class for fitting and interacting with Gaussian process models

#' @title Gaussian process models
#' 
#' @description Fitting (latent) Gaussian process models using a range of 
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
#' function to be used to fit the model. Currently only \code{gaussian}, 
#' \code{poisson} and \code{binomial} (only Bernoulli) are supported.
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
#' \code{gaussian(link = 'identity')} is full, direct inference, for
#' \code{binomial(link = 'logit')} and \code{binomial(link = 'probit')}
#' the default is full, Laplace inference (though note that only bernoulli)
#' data is handled at the moment.
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
#' # construct a poisson response variable
#' y2 <- rpois(n, exp(f))
#' 
#' # fit a GP model by Laplace approximation
#' # (note no observation error in this model)
#' m3 <- gp(y2, rbf('x'), df, poisson)
#' 
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



# user-facing predict function 
# add more options later to predicting with a new kernel,
# integrating over the point posterior, returning the full 
# posterior covariance etc.
#' @title predict method for Gaussian process models
#' 
#' @description Makes predictions from a fitted GP model in a similar way 
#' to \code{predict.glm} from GLMs
#' 
#' @param object a fitted GP model; an object of class \code{gp} produced by
#' \code{\link{gp}}. 
#'
#' @param newdata an optional data frame containing new data to predict to.
#' If \code{NULL}, the training data is used instead.
#' 
#' @param type the type of prediction required. The default (\code{'link'} 
#' predicts the value of the latent Gaussian process at the new locations,
#' \code{'response'} instead predicts on the scale of the response variable.
#' 
#' @param sd whether to report standard deviations of the prediction.
#' 
#' @param \dots further arguments to be passed to other methods, currently 
#' ignored
#' 
#' @return if \code{sd = FALSE} a vector of predictions. If 
#' \code{sd = TRUE}, a list with components: \code{fit} the 
#' predictions (as for \code{sd = FALSE}); \code{sd} estimated
#' standard deviations of the latent gaussians.
#' 
#' fitted gp object for which there aren't yet any associated 
#' functions. But there will be.
#' 
# @details
#' 
#' @export
#' @name predict.gp
#' 
#' @examples
#' # make some fake data
#' n <- 100  # observations
#' 
#' # dataframes
#' df <- data.frame(x = sort(runif(n, -5, 5)))
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
#' # make predictions from it (only the mean)
#' mu <- predict(m1, prediction_df)
#' 
#' plot(mu ~ prediction_df$x, type = 'l', ylim = range(y))
#' points(y ~ df$x, pch = 16)
#' 
#' # now predict with the standard deviations
#' res <- predict(m1, prediction_df, sd = TRUE)
#' 
#' # get upper and lower 95% credible intervals
#' upper <- res$fit + 1.96 * res$sd
#' lower <- res$fit - 1.96 * res$sd
#'
#' plot(mu ~ prediction_df$x, type = 'l',
#'  ylim = range(y, upper, lower))
#' lines(upper ~ prediction_df$x, lty = 2)
#' lines(lower ~ prediction_df$x, lty = 2)
#' points(y ~ df$x, pch = 16)



predict.gp <- function(object, newdata = NULL,
                       type = c('link', 'response'),
                       sd = FALSE,
                       ...) {
  
  # match the type of prediction
  type <- match.arg(type)
  
  # if new data isn't provided, predict back to the training data
  if (is.null(newdata)) {
    newdata <- object$data$training_data
  }
  
  # get the appropriate projector function
  projector <- getProjector(object$posterior)
  
  # project to the new data
  
  if (sd) {
    
    # if they want sds too, calculate the diagonal components
    res <- projector(object$posterior,
                     newdata,
                     variance = 'diag')
    
  } else {
    
    # otherwise don't calculate any variance components
    res <- projector(object$posterior,
                     newdata,
                     variance = 'none')
    
  }
  
  # N.B. could also add an option for the full posterior covariance matrix
  
  # calculate stuff
  
  if (type =='link') {
    
    # get the likelihood object
    # note this needs a y argument too, which can be NULL
    link <- object$likelihood$link
    
    # get answer on the link scale
    fit <- link(NULL, res$mu)
    
  } else {
    
    # otherwise, fit is on the link scale
    fit <- res$mu
    
  }
  
  # now work out what to return
  if (sd) {
    
    # if they want the sds
    ans <- list(fit = fit,
                sd = sqrt(res$var))
    
  } else {
    
    # otherwise, just the fit
    ans <- fit
    
  }
  
  # return the result
  return (ans)
  
}

# given a posterior object, return a compatible projector function
getProjector <- function (posterior) {
  
  # get inference name
  inference_name <- posterior$inference_name
  
  # get the corresponding projector name
  projector_name <- gsub('inference',
                         'project',
                         inference_name)
  
  # get the corresponding function
  projector <- get(projector_name)
  
  # return this
  return (projector)
  
}