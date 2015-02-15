getLikelihood <- function(likelihood, link) {
 
  # given strings giving the names of a likelihood and link,
  # return an object containing the derivative functions
  
  # get the name of the function to look for
  likelihood_string <- paste(likelihood,
                             link,
                             sep = '_')
  
  # check it exists
  if (exists(likelihood_string)) {
    
    # if so, get the master function
    evaluator <- get(likelihood_string)
    
  } else {
    
    # otherwise throw a nice error
    stop (paste0('Could not find a function to evaluate ',
                 likelihood,
                 ' with a ',
                 link,
                 ' link.'))
    
  }
  
  # define the functions for each derivative
  
  d0 <- function(y, f, ...) {
    evaluator(y, f, which = 'd0', ...)
  }
  
  d1 <- function(y, f, ...) {
    evaluator(y, f, which = 'd1', ...)
  }
  
  d2 <- function(y, f, ...) {
    evaluator(y, f, which = 'd2', ...)
  }
  
  # put it all in a list
  ans <- list(name = likelihood_string,
              d0 = d0,
              d1 = d1,
              d2 = d2)
  
  class(ans) <- 'likelihood'
  
  return (ans)
  
}