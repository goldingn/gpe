# likelihood_bernoulli_logit

likelihood_bernoulli_logit <- function(y, f, wt, which = c('d0', 'd1', 'd2', 'link'), ...) {
  # bernoulli log-likelihood (and its derivatives)
  # with the logit link function
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

    # hold tight, this is going to get messy...
    
    expfmy <- exp(f[pr] - y[pr])
    
    a <- expfmy / (expfmy + 1) ^ 2
    
    b <- (2 * exp(2 * f[pr] - 2 * y[pr])) / (expfmy + 1) ^ 3

    c <- exp(y[pr] - f[pr]) * (expfmy + 1) ^ 2
    
    ans[pr] <- (a - b) * c
    
    # for binary data...
    ans[npr] <- (y[npr] + 1) / 2 - plogis(f[npr])
    
  } else if (which == 'd2') {
    # second derivative
    
    # proportion data
 
    # convert to logistic scale
    y[pr] <- qlogis(y[pr])
    
    # even worse than the first derivative...
    expfmy <- exp(f[pr] - y[pr])
    exp2fmy <- exp(2 * (f[pr] - y[pr]))
    exp3fmy <- exp(3 * (f[pr] - y[pr]))
    expfmyp1 <- expfmy + 1
    expymf <- exp(y[pr] - f[pr])
    
    a <- (expfmy - (2 * exp2fmy) / expfmyp1) * expymf
    
    b <- (expfmy - (6 * exp2fmy) / expfmyp1 +
            (6 * exp3fmy) / expfmyp1 ^ 2) * expymf
    
    c <- (1 / expfmyp1) * (2 * expfmy - (4 * exp2fmy) / expfmyp1)
    
    ans[pr] <- -a + b + c
    
    # for binary data
    # get \pi as in Rasmussen & Williams for all non-proportion data
    pi_ <- plogis(f[npr])

    ans[npr] <- -pi_ * (1 - pi_)
    
  } else {
    
    # otherwise the link function *not* logged
    ans <- plogis(f)
    
  }
  
  # apply weights
  ans <- ans * wt
  
  return (ans)
  
}