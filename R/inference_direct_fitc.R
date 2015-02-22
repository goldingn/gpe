# inference_direct_fitc

# FITC inference for Gaussian likelihood and sparse GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# new_data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to evaluate the model
# inducing_data ia a dataframe giving the locations of inducing points
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
inference_direct_fitc <- function(y,
                                  data,
                                  kernel,
                                  likelihood,
                                  mean_function,
                                  inducing_data) {

  # NB likelihood is ignored
  
  # switch to ignoring new data,
  # creating a posterior object which can do predictions
  # later switch to efficient cholesky version
  
  # create a diagonal matrix with a small amount of noise for inducing points
  Ixx_noise <- diag(nrow(data)) * 10 ^ -6
  
  # evaluate prior mean function
  mn_pri_x <- mean_function(data)
  
  # get self kernels for old and new data
  # want self-variance on the old data (so one arg)
  Kxx <- kernel(data)
  
  # FITC components
  Kzz <- kernel(inducing_data)
  Kzzi <- solve(Kzz)
  Kxz <- kernel(data, inducing_data)
  Kzx <- t(Kxz)
  
  # (from Quinonero-candela & Rasmussen)
  Qxx <- Kxz %*% Kzzi %*% Kzx
  
  # calculate $\Lambda$
  Lambda <- diag(diag(Kxx - Qxx + Ixx_noise))
  # invert $\Lambda$ (it's diagonal, so just that bit)
  iLambda <- diag(1 / diag(Lambda))
  
  # calculate $\Sigma$
  Sigma <- solve(Kzz + Kzx %*% iLambda %*% Kxz)
  
  # log marginal likelihood
  # look into speeding up determinant and inversion
  # I think it's in Vanhatalo et al. sparse ... disease ...
  lZ <- - (log(det(Qxx + Lambda))) / 2 + 
    (t(y - mn_prior) %*% solve(Qxx + Lambda) %*% (y - mn_prior)) / 2 -
    (nrow(data) * log(2 * pi)) / 2
  
  # return posterior object
  posterior <- createPosterior(inference_name = 'inference_direct_exact',
                               lZ = lZ,
                               data = data,
                               kernel = kernel,
                               likelihood = likelihood,
                               mean_function = mean_function,
                               inducing_data = inducing_data,
                               Kzzi = Kzzi,
                               Kzx = Kzx,
                               Sigma = Sigma,
                               iLambda = iLambda,
                               y = y,
                               mn_prior = mn_prior)
  
  return (posterior)
  
}

# projection
project_direct_fitc <- function(posterior, new_data) {

  # prior mean over the test locations
  mn_prior_xp <- posterior$mean_function(new_data)
  
  # prior covariance over the test locations
  Kxpxp <- posterior$kernel(new_data)

  # FITC components
  Kxpz <- posterior$kernel(new_data,
                           posterior$inducing_data)
  Kzxp <- t(Kxpz)
  Qxpxp <- Kxpz %*% posterior$components$Kzzi %*% Kzxp
  
  mu <- Kxpz %*% posterior$components$Sigma %*%
    posterior$components$Kzx %*%
    posterior$components$iLambda %*%
    (posterior$components$y - posterior$components$mn_prior) +
    mn_prior_xp
  
  K <- Kxpxp - Qxpxp + Kxpz %*% posterior$components$Sigma %*% Kzxp
  
  # return both
  ans <- list(mu = mu,
              K = K)
  
  return (ans)
  
}
