# likelihood_bernoulli_logit

bernoulli_logit <- function(y, f, which = c('d0', 'd1', 'd2'), ...) {
  # bernoulli log-likelihood (and its derivatives)
  # with the logit link function
  # y is the observed data, either binary (0, 1) or proportion
  # f is the value of the correspnding latent gaussians
  # which determines whether to return the log-likelihood (d0)
  # or it's first or second derivative
  # \dots doesn't do anything, but is there for compatibility with other
  # likelihood which might use it
  
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
    
    # for proportion data, convert to univariate logistic
    # and get log-density
    y[pr] <- qlogis(y[pr])
    ans[pr] <- dlogis(y[pr],
                     f[pr],
                     log = TRUE)
    
    # for true binary data, get cumulative density of the unit gaussian
    # evaluated at f
    ans[npr] <- plogis(y[npr] * f[npr],
                      log.p = TRUE)
    
  } else if (which == 'd1') {
    # first derivative
    
    # for proportion data, convert to univariate Gaussian and evaluate
    y[pr] <- qlogis(y[pr])
    if (any(pr)) stop ('No derivatives for proportions with logit link. Yet.')
#    ans[pr] <- y[pr] - f[pr]
    
# for binary data...
    ans[npr] <- (y[npr] + 1) / 2 - plogis(f[npr])
    
  } else {
    # second derivative
    
    # all proportions get -1, apparently
    if (any(pr)) stop ('No derivatives for proportions with logit link. Yet.')
    #     ans[pr] <- -1
    
    # get \pi as in Rasmussen & Williams for all non-proportion data
    pi_ <- plogis(f[npr])

    ans[npr] <- -pi_ * (1 - pi_)
    
  }
  
  return (ans)
  
}