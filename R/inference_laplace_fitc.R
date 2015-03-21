# inference_laplace_fitc

# Inference using the laplace approximation and a sparse GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# inducing_data ia a dataframe giving the locations of inducing points
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
inference_laplace_fitc <- function(y,
                                   data,
                                   kernel,
                                   likelihood,
                                   mean_function,
                                   inducing_data,
                                   weights,
                                   verbose = verbose) {
  
  # apply mean function to get prior mean at observation locations
  mn_prior <- mean_function(data)
  
  # control parameters
  tol <- 10 ^ -12
  itmax = 50

  # number of observations
  n <- nrow(data)
  m <- nrow(inducing_data)
  
  # self kernel (with observation noise)
  Kxx <- kernel(data, data)
  
  # inducing data kernel (without observation noise)
  Kzz <- kernel(inducing_data)
  
  # projection kernels
  Kxz <- kernel(data, inducing_data)
  Kzx <- t(Kxz)
  
  # decomposition of the noiseless inducing points
  Lzz <- jitchol(Kzz)
  
  # diagonal elements of Qxx
  B <- forwardsolve(t(Lzz), Kzx)
  Qxx_diag <- colSums(B ^ 2)
  
  # diagonal elements of the full covariance function
  Kxx_diag <- diag(kernel(data, diag = TRUE))
  
  # diagonal elements of Lambda matrix
  Lambda_diag <- Kxx_diag - Qxx_diag
  
  # initialise a
  a <- rep(0, n)
  
  # set f to the prior
  f <- mn_prior
  
  # initialise loop
  obj.old <- Inf
  obj <- -sum(likelihood$d0(y, f, weights))
  it <- 0
  
  # start newton iterations
  while ((obj.old - obj) > tol & it < itmax) {
    
    # increment iterator and update objective
    it <- it + 1
    obj.old <- obj
    
    # get the negative log Hessian and its root
    W <- -(likelihood$d2(y, f, weights))
    rW <- sqrt(W)
    
    # difference between posterior mode and prior
    cf <- f - mn_prior
    
    # get components of the lambda
    Lah <- 1 + rW * Lambda_diag * rW
    rWKxz <- matrix(rep(rW, m), ncol = m) * Kxz
    mat1 <- matrix(rep(Lah, m), ncol = m)
    Sigma <- Kzz + t(rWKxz) %*% (rWKxz / mat1)
    # Sigma is A

    # enforce symmetry
    Sigma <- (Sigma + t(Sigma)) / 2
    
    # get direction of the posterior mode
    b <- W * cf + likelihood$d1(y, f, weights)
    Lb <- t(forwardsolve(t(jitchol(Sigma)), t(rWKxz / mat1)))
    mat2 <- rW * (Lambda_diag * b + t(B) %*% (B %*% b))
    adiff <- b - rW * (mat2 / Lah - Lb %*% (t(Lb) %*% mat2)) - a
    
    # make sure it's a vector, not a matrix
    dim(adiff) <- NULL
    
    # find optimum step size toward the mode using Brent's method
    res <- optimise(laplace_psiline_fitc,
                    interval = c(0, 2),
                    adiff = adiff,
                    a = a,
                    Lambda_diag = Lambda_diag,
                    B = B,
                    y = y,
                    d0 = likelihood$d0,
                    mn = mn_prior,
                    wt = weights)
    
    # move to the new posterior mode
    a <- a + res$minimum * adiff
    f <- Lambda_diag * a + t(B) %*% (B %*% a) + mn_prior

    # again, make it a vector
    dim(f) <- NULL
    
    # get the value of the objective
    obj <- laplace_psi(a = a,
                       f = f,
                       mn = mn_prior,
                       y = y,
                       d0 = likelihood$d0,
                       wt = weights)
    
  }
  
  # recompute hessian at mode and it's root
  W <- -(likelihood$d2(y, f, weights))
  rW <- sqrt(W)
  
  # rerun matrix algebra
  Lah <- 1 + rW * Lambda_diag * rW
  rWKxz <- matrix(rep(rW, m), ncol = m) * Kxz
  mat1 <- matrix(rep(Lah, m), ncol = m)
  Sigma <- Kzz + t(rWKxz) %*% (rWKxz / mat1)
  
  # enforce symmetry
  Sigma <- (Sigma + t(Sigma)) / 2

  # return marginal negative log-likelihood
  lZ <- sum(likelihood$d0(y, f, weights)) +
    0.5 * sum(log(Lah)) +
    -sum(log(diag(Lzz))) +
    sum(log(diag(jitchol(Sigma))))
  
  # get whitened version of inducing point function values
  f_z <- Kzx %*% a
  a_z <- backsolve(Lzz, forwardsolve(t(Lzz), f_z))[, 1]
  
  # return posterior object
  posterior <- createPosterior(inference_name = 'inference_laplace_fitc',
                               lZ = lZ,
                               data = data,
                               kernel = kernel,
                               likelihood = likelihood,
                               mean_function = mean_function,
                               inducing_data = inducing_data,
                               weights = weights,
                               Lzz = Lzz,
                               a_z = a_z,
                               Kzz = Kzz,
                               Kzx = Kzx,
                               Lambda_diag = Lambda_diag,
                               Sigma = Sigma,
                               f = f,
                               mn_prior = mn_prior,
                               W = W)
  
  # return a posterior object
  return (posterior)
  
}

# projection for fitc inference
# projection for FITC models
project_laplace_fitc <- function(posterior,
                                new_data,
                                variance = c('none', 'diag', 'matrix')) {
  
  # get the required variance argument
  variance <- match.arg(variance)
  
  # prior mean over the test locations
  mn_prior_xp <- posterior$mean_function(new_data)
  
  # prior covariance over the test locations
  Kxpxp <- posterior$kernel(new_data)
  
  # FITC components
  Kxpz <- posterior$kernel(new_data,
                           posterior$inducing_data)
  
  mu <- (Kxpz %*% posterior$components$a_z + mn_prior_xp)[, 1]
  
  if (variance == 'none') {
    
    # if mean only
    var <- NULL
    
  } else {
    
    # compute common variance components
    Kzxp <- t(Kxpz)
    Qxpxp <- Kxpz %*% solve(posterior$components$Kzz) %*% Kzxp
    K <- Kxpxp - Qxpxp + Kxpz %*% posterior$components$Sigma %*% Kzxp    
    
    if (variance == 'diag') {
      
      # if diagonal (elementwise) variance only
      
      # make this more efficient!
      var <- diag(K)
      
    } else {
      
      # if full variance
      var <- K
      
    }
  }
  
  # return both
  ans <- list(mu = mu,
              var = var)
  
  return (ans)
  
}
