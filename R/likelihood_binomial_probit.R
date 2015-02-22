# likelihood_binomial_probit

likelihood_binomial_probit <- function(y, f, which = c('d0', 'd1', 'd2'), ...) {
  # binomial log-likelihood (and its derivatives)
  # with the probit link function
  # y is the observed data, either binary (0, 1) or proportion
  # f is the value of the correspnding latent gaussians
  # which determines whether to return the log-likelihood (d0)
  # or it's first or second derivative
  # \dots doesn't do anything, but is there for compatibility with other
  # likelihood which might use it
  
  # check whether it's passed as two columns (binomial) or as one (bernoulli)
  if (ncol(y) == 2) {
    # if they want a proper binomial, throw an error as this isn't implemented
    # yet
    
    stop ('The binomial distribution with sample sizes greater than one is not \
yet implemented. Sorry about that. The Bernoulli distribution \
          is available however.')

    } else if (ncol(y) == 1) {
    
    # call the bernoulli likelihood
    ans <- likelihood_bernoulli_probit(y,
                                       f,
                                       which = which,
                                       ...)
    
    return (ans)
    
  } else {

    # if there's any other number, throw an error
    stop ('For the binomial distribution the response variable must have either one or two columns')
    
  }
}