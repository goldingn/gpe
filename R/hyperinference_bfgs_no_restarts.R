# hyperinference_bfgs_no_restarts

# gradient-free BFGS optimisation with no restarts
hyperinference_bfgs_no_restarts <- function(model) {
  
  # run BFGS optimisation on the model only once
  model <- optimizeModelBFGS(model, 1)
  
  return (model)
  
}