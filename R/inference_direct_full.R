# inference_direct_full

# Direct inference for Gaussian likelihood and full GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
inference_direct_full <- function(y,
                                   data,
                                   kernel,
                                   likelihood,
                                   mean_function,
                                   inducing_data,
                                   verbose = verbose) {
  
  # NB likelihood and inducing_data are ignored
  
  # apply mean function to get prior mean at observation locations
  mn_prior <- mean_function(data)
  
  # self kernel (with observation noise)
  K <- kernel(data, data)
  
  # its cholesky decomposition
  L <- jitchol(K, verbose = verbose)
  
  # uncorrelated latent vector a
  a <- backsolve(L, forwardsolve(t(L), (y - mn_prior)))

  # log marginal likelihood
  lZ <- - (t(y - mn_prior) %*% a)[1, 1] / 2 -
    sum(log(diag(L))) -
    (n * log(2 * pi)) / 2
  
  # return posterior object
  posterior <- createPosterior(inference_name = 'inference_direct_full',
                               lZ = lZ,
                               data = data,
                               kernel = kernel,
                               likelihood = likelihood,
                               mean_function = mean_function,
                               inducing_data = inducing_data,
                               mn_prior = mn_prior,
                               L = L,
                               a = a)
  
  # return a posterior object
  return (posterior)
  
}

# projection for full inference
project_direct_full <- function(posterior,
                                new_data) {
  
  
  # prior mean over the test locations
  mn_prior_xp <- posterior$mean_function(new_data)
  
  # projection matrix
  Kxxp <- posterior$kernel(posterior$data,
                           new_data)
  
  # its transpose
  Kxpx <- t(Kxxp)
  
  # test data self-matrix
  Kxpxp <- posterior$kernel(new_data)
  
  # get posterior mean
  mu <- Kxpx %*% posterior$components$a + mn_prior_xp
  
  # get posterior covariance
  v <- backsolve(posterior$components$L,
                 Kxxp,
                 transpose = TRUE)
  
  K <- Kxpxp - crossprod(v)
  
  # NB can easily modify this to return only the diagonal elements
  # (variances) with kernel(..., diag = TRUE)
  # calculation of the diagonal of t(v) %*% v is also easy:
  # (colSums(v ^ 2))
  
  # return both
  ans <- list(mu = mu,
              K = K)
  
  return (ans)
  
}