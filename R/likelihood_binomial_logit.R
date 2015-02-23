# likelihood_binomial_logit

likelihood_binomial_logit <- function(y, f, which = c('d0', 'd1', 'd2', 'link'), ...) {
  # binomial log-likelihood (and its derivatives)
  # with the logit link function
  # y is the observed data, either binary (0, 1) or proportion
  # f is the value of the correspnding latent gaussians
  # which determines whether to return the log-likelihood (d0)
  # or it's first or second derivative
  # \dots doesn't do anything, but is there for compatibility with other
  # likelihood which might use it
  
  # handle NULLs passed for the link function
  if (is.null(y)) {
    if (which == 'link') {
      
      ans <- likelihood_bernoulli_logit(y,
                                        f,
                                        which = which,
                                        ...)
    } else {
      stop ('y cannot be NULL unless the link function is being used')
    }
  }
  
  # check whether y is passed as a matrix
  if (is.matrix(y)) {
    
    if (ncol(y) == 2) {
      # if they want a proper binomial, throw an error as this isn't implemented
      # yet   
      
      stop ('The binomial distribution with sample sizes greater than one is not \
            yet implemented. Sorry about that. The Bernoulli distribution \
            is available however.')
      
    } else if (ncol(y) == 1) {
      
      # if it's a single variable as a matrix, convert to a vector
      y <- y[, 1]
      
    } else {
      
      # if there's any other number, throw an error
      stop ('For the binomial distribution the response variable must have either one or two columns')
      
    }
  }
  
  if (is.vector(y)) {
    # if it's a vector call the bernoulli likelihood
    ans <- likelihood_bernoulli_logit(y,
                                      f,
                                      which = which,
                                      ...)
  }
  
  return (ans)
  
}