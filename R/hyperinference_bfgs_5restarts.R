# hyperinference_bfgs_5restarts

# gradient-free BFGS optimisation with 5 restarts
hyperinference_bfgs_5_restarts <- function(model) {

  # run BFGS optimisation on the model 5 times with random restarts
  model <- optimizeModelBFGS(model, 5)
  
  return (model)

}