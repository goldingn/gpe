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
  
  link <- function(y, f, ...) {
    evaluator(y, f, which = 'link', ...)
  }
  
  # put it all in a list
  ans <- list(name = likelihood_string,
              d0 = d0,
              d1 = d1,
              d2 = d2,
              link = link)
  
  class(ans) <- 'likelihood'
  
  return (ans)
  
}

# sanity checks for data types
checkNonNegative <- function (response) {
  # throw an error if the response data provided is not a positive integer
  if (any(response < 0)) {
    stop("negative values are not allowed for this family")    
  }
}

getFamily <- function (likelihood) {
  # given a likelihood object, get the family object from which it came
  
  # get the name
  likelihood_name <- likelihood$name
  
  # split it up
  likelihood_split <- strsplit(likelihood_name, '_')[[1]]
  
  # split out the family and link names
  family_name <- likelihood_split[2]
  link_name <- likelihood_split[3]
  
  # get the family object
  family <- eval(parse(text = paste0(family_name,
                                       '("',
                                       link_name,
                                       '")')))
  
  # return it
  return (family)
  
}

checkNonNegative <- function (response) {
  # throw an error if the response data provided is not a positive integer
  if (any(response < 0)) {
    stop ("negative values not allowed in the response variable with this family")    
  }
}

checkUnitInterval <- function (response) {
  # throw an error if the response data is not on the unit interval
  if (any(response < 0 | response > 1)) {
    stop ("response values greater than one or less than zero are not allowed for this family")    
  }

  # issue a warning if any of it is non-integer(ish)
  if (any(response > 0 & response < 1)) {
    warning ("non-integer values in the response variable, these will be treated as having an infinite sample size")    
  }
  
}

