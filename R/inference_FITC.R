# FITC inference for Gaussian likelihood and sparse GP

# y is response data
# data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to train the model
# new_data is a dataframe containing columns on which the kernel acts,
#   giving observations on which to evaluate the model
# inducing_data ia a dataframe giving the locations of inducing points
# kernel is a gpe kernel object
# mean_function is a gpe mean function (for now just an R function)
infFITC <- function(y,
                     data,
                     new_data,  # get rid of this!
                     inducing_data,
                     kernel,
                     mean_function) {
  
  # switch to ignoring new data,
  # creating a posterior object which can do predictions
  
  # evaluate prior mean function
  mn_pri_x <- meanfunction(data)
  mn_pri_xp <- meanfunction(newdata)

  # get self kernels for old and new data
  # want self-variance on the old data (so one arg)
  Kxx <- kernel(data)
  # but not on the new data as we want expected not observed
  # (so both args)
  Kxpxp <- kernel(newdata, newdata)
  
  # FITC components
  Kzz <- kernel(inducingdata)
  Kzzi <- solve(Kzz)
  Kxz <- kernel(data, inducingdata)
  Kzx <- t(Kxz)
  Kxpz <- kernel(newdata, inducingdata)
  Kzxp <- t(Kxpz)
  
  # (from Quinonero-candela & Rasmussen)
  Qxx <- Kxz %*% Kzzi %*% Kzx
  Qxpxp <- Kxpz %*% Kzzi %*% Kzxp

  # calculate $\Lambda$
  Lambda <- diag(diag(Kxx - Qxx))
  # invert $\Lambda$ (it's diagonal, so just that bit)
  iLambda <- diag(1 / diag(Lambda))
  
  # calculate $\Sigma$
  Sigma <- solve(Kzz + Kzx %*% iLambda %*% Kxz)
  
  mn <- Kxpz %*% Sigma %*% Kzx %*% iLambda %*% (y - mn_pri_x) +
    mn_pri_xp
  K <- Kxpxp - Qxpxp + Kxpz %*% Sigma %*% Kzxp
  
  ans <- list(mn = mn, K = K)
  
  return (ans)

}
