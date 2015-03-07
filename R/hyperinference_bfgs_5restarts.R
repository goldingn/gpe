# hyperinference_bfgs_5restarts

# gradient-free BFGS optimisation with 5 restarts
hyperinference_bfgs_5_restarts <- function(model) {

  # run BFGS optimisation on the model 5 times with random restarts
  # each set of parameters sampled from a random normal with mean
  # given by the current parameter set and a standard deviation of 10
  model <- optimizeModelBFGS(model,
                             restarts = 5,
                             sampling_sd = 10)
  
  return (model)

}