# miscellaneous utility functions for all classes

clone <- function(x) {
  # clone an object, making sure the new object has a different environment
  
  # get the environment of x
  e_old <- environment(x)
  
  # turn into as a list
  e_old_list <- as.list(e_old, all.names = TRUE)
  
  # create a new environment
  e_new <- as.environment(e_old_list)
  
  # copy x
  environment(x) <- e_new
  
  return (x)
  
}