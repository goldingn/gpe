# likelihood_bernoulli_probit

likelihood_bernoulli_probit <- function(y, f, wt, which = c('d0', 'd1', 'd2', 'link'), ...) {
  # bernoulli log-likelihood (and its derivatives)
  # with the probit link function
  # y is the observed data, either binary (0, 1) or proportion
  # f is the value of the correspnding latent gaussians
  # which determines whether to return the log-likelihood (d0)
  # or it's first or second derivative
  # \dots doesn't do anything, but is there for compatibility with other
  # likelihood which might use it

  # check response variable
  checkUnitInterval(y)
  
  # check which derivative is needed
  which <- match.arg(which)
  
  # repeat y if needed
  if(length(y) != length(f)) y <- rep(y, length(f))
  
  # find observations which are actually probabilities
  # (can be used to handle proportion data with unknown sample size)

  # proportions
  pr <- y > 0 & y < 1
  
  # not proportions
  npr <- !pr

  # switch true binary data from 0, 1 to -1, 1
  y[npr] <- 2 * y[npr] - 1
  
  # create an empty vector for the results
  ans <- vector('numeric', length(y))
  
  # straight likelihood case
  if (which == 'd0') {
  
    # for proportion data, convert to univariate Gaussian
    # and get log-density
    y[pr] <- qnorm(y[pr])
    ans[pr] <- dnorm(y[pr],
                     f[pr],
                     log = TRUE)
    
    # for true binary data, get cumulative density of the unit gaussian
    # evaluated at f
    ans[npr] <- pnorm(y[npr] * f[npr],
                      log.p = TRUE)
  
  } else if (which == 'd1') {
    # first derivative
    
    # for proportion data, convert to univariate Gaussian and evaluate
    y[pr] <- qnorm(y[pr])
    ans[pr] <- y[pr] - f[pr]
    
    # for binary data...
    ans[npr] <- y[npr] * dnorm(f[npr]) / pnorm(y[npr] * f[npr])
    
  } else if (which == 'd2') {
    # second derivative
    
    ans[pr] <- -1
    
    # the binary data gets something more complex
    a <- dnorm(f[npr]) ^ 2 / pnorm(y[npr] * f[npr]) ^ 2
    b <- y[npr] * f[npr] * dnorm(f[npr]) / pnorm(y[npr] * f[npr])
    ans[npr] <- -a - b

  } else {
    
    # otherwise the link, *not* on the log scale
    ans <- pnorm(f)
    
  }
  
  # apply weights
  ans <- ans * wt
  
  return (ans)

}