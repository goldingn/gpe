# hyperinference_utility

# utility functions for inference over the hyperparameters of a gp model
# (the kernel parameters)

nlmlObjective <- function(pars, model) {
  # return the negative log marginal likelihood of a GP model with kernel 
  # parameters pars
  
  # get the model kernel
  kernel <- model$kernel
  
  # update the parameters
  kernel <- gpe:::setParVec(kernel, par_vec = pars)
  
  # fit the model with the new parameters
  model_new <- gpe:::updateModel(model = model, kernel = kernel)
  
  # extract the negative log marginal likelihood
  nlml <- -model_new$posterior$lZ
  
  # return this
  return (nlml)
  
}


optimizeModelBFGS <- function (model, restarts = 1) {
  # optimize the hyperparameters of a model by gradient-free BFGS,
  # optionally with random restarts.
  # Note that the non-convexity of the likelihood means that running this without
  # random restarts is probably a terrible idea (at least without a strong 
  # hyperprior), particuarly on Gaussian-response data.

  # get current kernel
  kernel <- model$kernel
  
  # get current parameters of that kernel
  pars <- getParVec(kernel)
  
  # get a list of random parameters for each restart
  pars_list <- replicate(restarts,
                         pars + rnorm(length(pars), 0, 10),
                         simplify = FALSE)
  
  # run the optimiser for each of these start values
  # sequrntially for the time being
  opt_list <- lapply(pars_list,
                     optim,
                     fn = nlmlObjective,
                     gr = NULL,
                     model = model)
  
  # get the objective functions of these
  objectives <- lapply(opt_list, function(x) x$value)
  
  # find the best model
  opt <- opt_list[[which.min(objectives)]]
  
  # reconstruct the new kernel
  kernel_new <- setParVec(kernel, opt$par)
  
  # fit the new model
  model_new <- updateModel(model, kernel = kernel_new)
  
  return (model_new)
  
}
