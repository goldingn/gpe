# likelihood_poisson_log

likelihood_poisson_log <- function(y, f, which = c('d0', 'd1', 'd2', 'link'), ...) {
  # poisson log-likelihood (and its derivatives)
  # with the log link function
  # y is the observed data, must be counts
  # f is the value of the corresponding latent gaussians
  # which determines whether to return the log-likelihood (d0)
  # or it's first or second derivative
  # \dots doesn't do anything, but is there for compatibility with other
  # likelihood which might use it
  
  # check which derivative is needed
  which <- match.arg(which)
  
  # repeat y if needed
  if(length(y) != length(f)) y <- rep(y, length(f))
  
  # create an empty vector for the results
  ans <- vector('numeric', length(y))
  
  # straight likelihood case
  if (which == 'd0') {

    # get log-density
    ans <- dpois(y,
                 exp(f),
                 log = TRUE)
    
  } else if (which == 'd1') {
    
    # first derivative
    ans <- y / exp(f) - 1
    
  } else if (which == 'd2') {

    # second derivative
    ans <- -y / (exp(f) ^ 2)
    
  } else {
    
    # otherwise link
    ans <- exp(f)
    
  }
  
  return (ans)
  
}