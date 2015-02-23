# inference_laplace_full

# Exact inference using the laplace approximation and a full GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
inference_laplace_full <- function(y,
                                   data,
                                   kernel,
                                   likelihood,
                                   mean_function,
                                   inducing_data,
                                   verbose = verbose) {
  
  # NB inducing_data is ignored
  
  # apply mean function to get prior mean at observation locations
  mn_prior <- mean_function(data)
  
  # control parameters
  tol <- 10 ^ -12
  itmax = 50
  
  # self kernel (with observation noise)
  Kxx <- kernel(data, data)
  
  # number of observations
  n <- nrow(data)

  # diagonal matrix
  eye <- diag(n)
  
  # initialise a
  a <- rep(0, n)
  
  # set f to the prior
  f <- mn_prior
  
  # initialise loop
  obj.old <- Inf
  obj <- -sum(likelihood$d0(y, f))
  it <- 0
  
  # start newton iterations
  while ((obj.old - obj) > tol & it < itmax) {

    # increment iterator and update objective
    it <- it + 1
    obj.old <- obj
    
    # get the negative log Hessian and its root
    W <- -(likelihood$d2(y, f))
    rW <- sqrt(W)
    
    # difference between posterior mode and prior
    cf <- f - mn_prior
    
    # get cholesky factorisation
    L <- jitchol(rW %*% t(rW) * Kxx + eye)
    
    # get direction of the posterior mode
    b <- W * cf + likelihood$d1(y, f)
    mat2 <- rW * (Kxx %*% b)
    adiff <- b - rW * backsolve(L, forwardsolve(t(L), mat2)) - a 
    
    # make sure it's a vector, not a matrix
    dim(adiff) <- NULL
    
    # find optimum step size toward the mode using Brent's method
    res <- optimise(laplace_psiline_full,
                    c(0, 2),
                    adiff,
                    a,
                    Kxx,
                    y,
                    likelihood$d0,
                    mn_prior)
    
    # move to the new posterior mode
    a <- a + res$minimum * adiff
    f <- Kxx %*% a + mn_prior
    
    obj <- laplace_psi(a,
                       f,
                       mn_prior,
                       y,
                       likelihood$d0)
    
  }

  # recompute hessian at mode
  W <- -(likelihood$d2(y, f))
  
  # return marginal negative log-likelihood
  lZ <- -(a %*% (f - mn_prior))[1, 1] / 2 -
    sum(likelihood$d0(y, f)) +
    sum(log(diag(L)))
  
  # return posterior object
  posterior <- createPosterior(inference_name = 'inference_laplace_full',
                               lZ = lZ,
                               data = data,
                               kernel = kernel,
                               likelihood = likelihood,
                               mean_function = mean_function,
                               inducing_data = inducing_data,
                               mn_prior = mn_prior,
                               L = L,
                               a = a,
                               W = W)
  
  # return a posterior object
  return (posterior)
  
}

# projection for full inference
project_laplace_full <- function(posterior,
                                new_data,
                                variance = c('none', 'diag', 'matrix')) {
  
  # get the required variance argument
  variance <- match.arg(variance)
  
  # prior mean over the test locations
  mn_prior_xp <- posterior$mean_function(new_data)
  
  # projection matrix
  Kxxp <- posterior$kernel(posterior$data,
                           new_data)
  
  # its transpose
  Kxpx <- t(Kxxp)
  
  # get posterior mean
  mu <- Kxpx %*% posterior$components$a + mn_prior_xp
  
  # NB can easily modify this to return only the diagonal elements
  # (variances) with kernel(..., diag = TRUE)
  # calculation of the diagonal of t(v) %*% v is also easy:
  # (colSums(v ^ 2))
  
  if (variance == 'none') {
    
    # if mean only
    var <- NULL
    
  } else {
    
    # compute common variance components
    rW <- sqrt(as.vector(posterior$components$W))
    
    # get posterior covariance
    v <- backsolve(posterior$components$L,
                   rW * Kxxp,
                   transpose = TRUE)
    

    if (variance == 'diag') {
      
      # if diagonal (elementwise) variance only
      # diagonal matrix of the prior covariance on xp
      Kxpxp_diag <- posterior$kernel(new_data, diag = TRUE)
      
      # diagonal elements of t(v) %*% v
      vtv_diag <- colSums(v ^ 2)
      
      # diagonal elements of the posterior
      K_diag <- diag(Kxpxp_diag) - vtv_diag
      
      var <- K_diag
      
    } else {
      
      # if full variance
      
      # prior covariance on xp
      Kxpxp <- posterior$kernel(new_data)
      
      # posterior covariance on xp
      K <- Kxpxp - crossprod(v)

      var <- K
      
    }
  }
  
  # return both
  ans <- list(mu = mu,
              var = var)
  
  return (ans)
  
}