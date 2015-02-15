# Exact inference for Gaussian likelihood and full GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
infExact <- function(y,
                     data,
                     kernel,
                     mean_function) {
  
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
  
  # return list
  ans <- list(mu = mu,  # posterior mean
              K = K)  # posterior covariance matrix
  
  return (ans)

}
