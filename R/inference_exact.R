# exact inference for Gaussian likelihood and full GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# newdata is a dataframe containing columns on which the kernel acts,
#   giving observations on which to evaluate the model
# kernel is a gpe kernel object
# meanfunction is a gpe mean function (for now just an R function)
infExact <- function(y,
                     data,
                     newdata,
                     kernel,
                     meanfunction) {
  
  mn_prior <- meanfunction(data)
  mn_posterior <- meanfunction(newdata)
  Kxx <- kernel(data)
  Kxxi <- solve(Kxx)
  Kxxp <- kernel(data, newdata)
  Kxpx <- t(Kxxp)
  
  mn <- Kxpx %*% Kxxi %*% (y - mn_prior) + mn_posterior
  K <- Kxpx %*% Kxxi %*% Kxxp
  ans <- list(mn = mn, K = K)
  return (ans)
}
