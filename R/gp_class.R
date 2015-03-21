# functions for the gp class for fitting and interacting with Gaussian 
# process models

#' @name gp
#' @rdname gp
#'
#' @title Gaussian process models
#' 
#' @description Fitting, summarizing and predicting from Gaussian process 
#' models using a range of inference methods.
#' 
NULL

#' @rdname gp
#' @export
#' @param formula an object of class \code{\link{formula}} giving a symbolic 
#' description of the model to be fitted - the right hand side must be a 
#' valid kernel object or the commands to construct one.
#' 
#' @param data a data frame containing the covariates against which to model 
#' the response variable. This must have the same number of rows as 
#' \code{response} and contain named variables matching those referred to by
#' \code{kernel}.
#' 
#' @param family a \code{\link{family}} object giving the likelihood and link 
#' function to be used to fit the model. Currently only \code{gaussian}, 
#' \code{poisson} and \code{binomial} (only Bernoulli) are supported.
#' 
#' @param weights an optional vector of 'prior weights' to be used in 
#' the fitting process. Should be NULL or a numeric vector.
#' 
#' @param mean_function an optional function specifying the prior over the mean
#' of the gp, in other words a 'first guess' at what the true function is.
#' This must act on a dataframe with named variables matching some of those in 
#' \code{data} and return a vector giving a single value for each row in the 
#' dataframe. Note that this function must return a prediction on the scale of 
#' the link, rather than the response. If \code{NULL} then a prior mean is 
#' assumed to be 0 for all observations.
#' 
#' @param inducing_data an optional dataframe containing the locations of 
#' inducing points to be used when carrying out sparse inference (e.g. FITC).
#' This must contain variables with names matching those referenced by 
#' \code{kernel} and \code{mean_function}. This should have fewer rows than 
#' \code{data} and \code{response}. 
#' 
#' @param inference a string specifying the inference method to be used to 
#' estimate the values of the latent parameters. If \code{'default'} an 
#' appropriate method is picked for the likelihood specified. See details 
#' section for the list of default inference methods.
#' 
#' @param hyperinference the method to be used for inference on the 
#' hyperparameters (parameters of the kernel). \code{BFGS} carries out
#' straightforward optimisation starting from the current kernel 
#' parameters using gradient-free BFGS. Because the likelihood surface is
#' rarely convex, this is is not advised for general use. \code{BFGSrestarts}
#' runs gradient-free BFGS 5 times, each starting with a different randomly 
#' chosen set of parameters, which might be a bit better. \code{none} does
#' no inference on the hyperparameters. \code{default} currently switches to 
#' \code{none} in all cases, though in future it may depend on other arguments.
#' 
#' @param verbose whether to return non-critical information to the user 
#' during model fitting.
#' 
#' @return A fitted gp object for which there aren't yet any associated 
#' functions. But there will be.
#' 
#' @details The default inference method for a model with the family 
#' \code{gaussian(link = 'identity')} is full direct inference (\code{'full'}),
#' for \code{binomial(link = 'logit')} and \code{binomial(link = 'probit')}
#' the default is full Laplace inference (\code{'Laplace'}; though note that 
#' only Bernoulli data is handled at the moment). Sparse inference can be 
#' carried out by specifying \code{inference = 'FITC'}, this is currently only
#' available for a model with a Gaussian likelihood.
#' 
#' @examples
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
#' # fit a full (non-sparse) GP model (without updating the hyperparameters) 
#' # as this is the default. Notice we add the observation error to the kernel.
#' m1 <- gp(y ~ rbf('x') + iid(), df, gaussian)
#' 
#' # fit another with FITC sparsity
#' m2 <- gp(y ~ rbf('x') + iid(), df, gaussian, inference = 'FITC', 
#'          inducing_data = df)
#'          
#' # summary stats, other associated functions still to come
#' 
#' # construct a poisson response variable
#' y2 <- rpois(n, exp(f))
#' 
#' # fit a GP model by Laplace approximation
#' # (note no observation error in this model)
#' m3 <- gp(y2 ~ rbf('x'), df, poisson)
#' 
gp <- function(formula,
               data,
               family = gaussian,
               weights = NULL,
               mean_function = NULL,
               inducing_data = NULL,
               inference = c('default', 'full', 'FITC', 'Laplace', 'LaplaceFITC'),
               hyperinference = c('default', 'BFGS', 'BFGSrestarts', 'none'),
               verbose = FALSE) {
  
  # catch the call
  call <- match.call()
  
  # get response and kernel from the formula
  response <- getResponse(formula, data)
  kernel <- getKernel(formula)
  
  # get or check the weights
  weights <- getWeights(weights, length(response))
  
  # get or check the mean function
  mean_function <- getMeanFunction(mean_function, data)
  
  # get the likelihood
  likelihood <- getLikelihood(family)
  
  # match the inference option
  inference <- match.arg(inference)
  
  # get inference function
  inference <- getInference(inference,
                            likelihood,
                            inducing_data)
  
  # match the hyperparameter inference
  hyperinference_name <- match.arg(hyperinference)
  
  # get hyperparameter inference function
  hyperinference <- getHyperinference(hyperinference_name)
  
  # check it's a valid kernel
  checkKernel(kernel)  
  
  # do inference
  posterior <- inference(response,
                         data,
                         kernel,
                         likelihood,
                         mean_function,
                         inducing_data,
                         weights,
                         verbose = verbose)
  
  model <- list(call = call,
                likelihood = likelihood,
                kernel = kernel,
                mean_function = mean_function,
                data = list(response = response,
                            training_data = data,
                            inducing_data = inducing_data),
                posterior = posterior,
                hyperinference = hyperinference_name,
                verbose = verbose)
  
  class(model) <- 'gp'
  
  # do hyperparameter inference
  model <- hyperinference(model)
  
  return (model)
  
}

#' @rdname gp
#' 
#' @export
#' 
#' @param x an object of class \code{gp}, constructed by the function \code{gp}
#' giving a fitted gaussian process model
#' 
#' @param \dots additional arguments for compatibility with generic functions
#'  
#' @examples
#' 
#' print(m3)
#' m3
#' 
print.gp <- function (x, ...) {
  # basic print function for a gp object (fitted gp model)
  
  cat('A Gaussian process model fitted against',
      nrow(x$data$training_data),
      'observations\n')
  
  cat('\tcall: \t\t\t',
      deparse(x$call),
      '\n')
  
  cat('\tkernel:',
      capture.output(print(x$kernel)),
      '\n')
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
#' # fit a full (non-sparse) GP model (without updating the hyperparameters)
#' # with rbf kernel and observation error
#' m1 <- gp(y ~ rbf('x') + iid(), df, gaussian)
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
  
  if (type =='response') {
    
    # get the likelihood object
    # note this needs a y argument too, which can be NULL
    link <- object$likelihood$link
    
    # get answer on the link scale
    fit <- link(res$mu)
    
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

is.gp <- function (x) {
  
  # test whether x is a gp object
  ans <- inherits(x, "gp")
  
  # return the answer
  return (ans)
  
}


getResponse <- function(formula, data) {
  # given a formula and a dataframe, extract the response variable
  # first looking in the dataframe, then the parent environment
  
  # get name of reponse variable
  response_name <- as.character(formula[[2]])
  
  # names of entries in dataframe
  data_names <- names(data)
  
  # if it's in there...
  if (response_name %in% data_names) {
    
    # get it
    response <- data[, match(response_name, data_names)]
    
  } else {
    
    # otherwise look up for it
    response <- get(response_name, envir = environment(formula))
    
  }
  
  # return
  return (response)
  
}

getKernel <- function (formula) {
  # given a formula, fetch the kernel object
  
  # get the string from the formula
  kernel_string <- deparse(formula[[3]])
  
  # first try evaluating it
  kernel <- eval(parse(text = kernel_string))
  
  # if it's not a valid kernel
  if (!is.kernel(kernel)) {
    
    # fetch the object with that name from the environment of the formula
    kernel <- get(kernel_string, envir = environment(formula))
    
  }
  
  # check it's a kernel object
  checkKernel(kernel)
  
  # and return
  return (kernel)
  
}

checkWeights <- function (weights, n) {
  # check the weights argument and throw a nice error if invalid
  
  # must be a vector
  if (!is.vector(weights)) {
    stop ('weights must be a vector, or NULL')
  }
  
  # must have length n
  if (length(weights) != n) {
    stop ('weights must have the same length as the response variable')
  }
  
  # must be numeric
  if (!is.numeric(weights)) {
    stop ('weights must be numeric')
  }
  
  # must be finite
  if (any(!is.finite(weights))) {
    stop ('weights must be finite')
  }
  
  # can't be negative
  if (any(weights < 0)) {
    stop ('negative weights not allowed')
  }
  
}


getWeights <- function (weights, n) {
  # check the weights argument and either return valid weights
  # or throw a nice error
  
  if (is.null(weights)) {
    
    # if null, return ones
    weights <- rep(1, n)
    
  } else {
    
    # otherwise check the weights passed are valid
    checkWeights(weights, n)
    
  }
  
  # return the weights
  return (weights)
  
}

checkMeanFunction <- function(mean_function, data) {
  # check that the specified mean function is valid
  
  # get the number of elements in the dataframe
  n <- nrow(data)
  
  # evaluate it on the dataframe
  mn <- mean_function(data)
  
  # mn must be finite
  if (any(!is.finite(mn))) {
    stop ('mean_function returned non-finite values')
  }
  
  # must be a vector
  if (!is.vector(mn)) {
    stop ('mean_function must return a vector')
  }
  
  # must have length n
  if (length(mn) != n) {
    stop ('mean_function must return a vector of the same length as
          the dataframe against which it is evaluated')
  }
  
  # must be numeric
  if (!is.numeric(mn)) {
    stop ('mean_function must return numeric values')
  }
  
  }

getMeanFunction <- function (mean_function, data) {
  # check the weights argument and either return valid weights
  # or throw a nice error
  
  if (is.null(mean_function)) {
    
    # if null, return ones
    mean_function <- function (data) {
      rep(0, nrow(data))
    }
    
  } else {
    
    # otherwise check the weights passed are valid
    checkMeanFunction(mean_function, data)
    
  }
  
  # return the weights
  return (mean_function)
  
}
