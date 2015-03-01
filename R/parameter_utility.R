# parameter_utility

# utility functions for parameter class objects

# boilerplate function to aid construction of parameter objects
# type is the parameter type (support), value is the initial value it takes
# true2cont is a function moving from value to a continuous transformation
# cont2true is the inverse function
# test is a function which returns TRUE if the function falls within the 
# correct constraints
createParameterConstructor <- function (type,
                                        value,
                                        true2cont,
                                        cont2true,
                                        test) {
  
  # create model object and initialize
  object <- list(type = type,
                 value = value,
                 true2cont = true2cont,
                 cont2true = cont2true,
                 test = test)
  
  # create function to return
  # it either evaluates as the true value (default)
  # or its continuous transform
  ans <- function(continuous = FALSE) {
    if (continuous) {
      object$true2cont(object$value)
    } else {
      object$value
    }
  }
  
  class(ans) <- 'parameter'
  
  return (ans)
}

# get the cont2true function of a parameter object
cont2true <- function (parameter) {
  
  # check it's valid
  checkParameter(parameter)
  
  # get the function
  ans <- environment(parameter)$object$cont2true
  
  # return
  return (ans)
  
}

# get the true2cont function of a parameter object
true2cont <- function (parameter) {
  
  # check it's valid
  checkParameter(parameter)
  
  # get the function
  ans <- environment(parameter)$object$true2cont
  
  # return
  return (ans)
  
}

# get the test function of a parameter object
test <- function (parameter) {
  
  # check it's valid
  checkParameter(parameter)
  
  # get the function
  ans <- environment(parameter)$object$test
  
  # return
  return (ans)
  
}

# check that x is a parameter
checkParameter <- function (x) {
  if (!is.parameter(x)) {
    stop (paste0('a parameter object is required, but an object of class ',
                 paste0(class(x), collapse = ', '),
                 ' was passed instead'))
  }
}


fun2string <- function (fn) {
  # return a string representing a function without ugly bits like
  # environments being displayed
  
  # copy over the function
  tmp_fn <- fn
  
  # assign it to the global environment (doesn't get reported on print)
  environment(tmp_fn) <- globalenv()
  
  # capture as a string
  fn_string <- capture.output(print(tmp_fn))
  
  # remove tmp_fn
  rm(tmp_fn)
  
  # return this
  return (fn_string)
  
}
