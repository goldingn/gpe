# inference_utility

# default mean function
zeroes <- function (data) rep(0, nrow(data))

# create a posterior object (doesn't actually do much ...)
createPosterior <- function (mu,
                             K,
                             X,
                             kernel,
                             inference) {
  
  # combine these things into a posterior
  ans <- list(mu = mu,
              K = K,
              X = X,
              kernel = kernel,
              inference = inference)
  
  # set it to a posterior class
  class(ans) <- c(inference, 'posterior')
  
  return (ans)
  
}
