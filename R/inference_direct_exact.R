# inference_direct_exact

# Exact inference for Gaussian likelihood and full GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
inference_direct_exact <- function(y,
                     data,
                     kernel,
                     mean_function,
                     inducing_data = NULL) {
  
  # mean function
  mn_prior <- mean_function(data)
  
  # self kernel and inverse
  Kxx <- kernel(data)
  Kxxi <- solve(Kxx)
  
  # predict to itself
  # (different from self kernel, due to observation variance)
  Kxxp <- kernel(data, data)
  Kxpx <- t(Kxxp)
  
  # get posterior mean and covariance matrix
  mu <- Kxpx %*% Kxxi %*% (y - mn_prior)
  K <- Kxpx %*% Kxxi %*% Kxxp
  
  # return posterior object
  posterior <- createPosterior(mu = mu,
                               K = K,
                               X = data,
                               kernel = kernel,
                               inference = 'inference_direct_exact')

  # need to add log marginal likelihood
  
  # return a posterior object
  return (posterior)

}

# projection
project_direct_exact <- function(posterior, new_data) {
  Kxpx <- posterior$kernel(new_data,
                           posterior$X)
  
  # project the posterior mean
  mu <- Kxpx %*% posterior$mu
  
  # get the posterior covariance
  Kxpxp <- posterior$kernel(new_data,
                            new_data)
  
  # return both
  ans <- list(mu = mu,
              K = Kxpxp)
  
  return (ans)

}
