# hyperinference_bfgs_no_restarts

# gradient-free BFGS optimisation with no restarts
hyperinference_bfgs_no_restarts <- function(model) {
  
  # run BFGS optimisation on the model only once
  model <- optimizeModelBFGS(model,
                             restarts = 1,
                             sampling_sd = 1)
  
  return (model)
  
}