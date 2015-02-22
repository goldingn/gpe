getLikelihood <- function(family) {

  # given the family argument passed to gp, parse and check the input and
  # return an object containing the derivative functions
  
  # adapted this from glm():
  if (is.character(family)) 
    
    # look up two environments (assume getLikelihood is directly within gp)
    # if that changes, change n here
    family <- get(family,
                  mode = "function",
                  envir = parent.frame(n = 2))

  if (is.function(family)) 
  
    family <- family()
  
  if (is.null(family$family)) {
  
    # does this ever get called?
    print(family)
    stop("'family' not recognized")
  
  } 
  
  # family should now be a family object
  # get the likelihood and link names
  
  likelihood <- family$family
  link <- family$link
  
  # get the name of the function to look for
  likelihood_string <- paste('likelihood',
                             likelihood,
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