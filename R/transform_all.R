# class of transformation constructors, similar to kernel mechanism

# boilerplate function to aid construction of transformation
# functions.
# type is the transform name, value is the initial value it takes
# true2cont is a function moving from value to a continuous transformation
# cont2true is the inverse function
# test is a function which returns TRUE if the function falls within the 
# correct constraints
createTransformConstructor <- function (type,
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
  ans <- function(cont = FALSE) {
    if (cont) {
      object$true2cont(object$value)
    } else {
      object$value
    }
  }
 
  class(ans) <- 'transform'

  return (ans)
}

# get the cont2true function of a transform object
cont2true <- function (transform) {
  
  # check it's valid
  checkTransform(transform)
  
  # get the function
  ans <- environment(transform)$object$cont2true
  
  # return
  return (ans)
  
}

# get the true2cont function of a transform object
true2cont <- function (transform) {
  
  # check it's valid
  checkTransform(transform)
  
  # get the function
  ans <- environment(transform)$object$true2cont
  
  # return
  return (ans)
  
}

# get the test function of a transform object
test <- function (transform) {
  
  # check it's valid
  checkTransform(transform)
  
  # get the function
  ans <- environment(transform)$object$test
  
  # return
  return (ans)
  
}

# test whether x is a transform
is.transform <- function (x) {
  ans <- inherits(x, 'transform')
  return (ans)
}

# check that x is a transform
checkTransform <- function (x) {
  if (!is.transform(x)) {
    stop (paste0('a transform object is required, but an object of class ',
                 paste0(class(x), collapse = ', '),
                 ' was passed instead'))
  }
}

# function to update value of a transform
# transform is the transform,
# new_value is the new value
# cont specifyies whether new_value is provided on the 
update.transform <- function(transform, new_value, cont = FALSE) {
  
  # check it's a valid transform
  checkTransform(transform)
  
  # if the value is provided on the continuous scale
  if (cont) {

    # switch it to the true scale
    new_value <- cont2true(transform)(new_value)

  }
  
  # test the new value
  if (!test(transform)(new_value)) {
    stop ('new_value is not within the constraints of this parameter')
  }
  
  # update the value
  environment(transform)$object$value <- new_value
  
  # return the transform
  return (transform)
}

# print function for transform objects
print.transform <- function (transform) {
  
  # first check it's a valid transform object
  checkTransform(transform)
  
  # get type
  type <- environment(transform)$object$type
  
  # print all of these things
  cat('transform object of type:',
      type,
      '\n\ttest function:\t\t\t\t\t\t\t',
      fun2string(test(transform)),
      '\n\tcont2true function:\t\t',
      fun2string(cont2true(transform)),
      '\n\ttrue2cont function:\t\t',
      fun2string(true2cont(transform)),
      '\n')
  
}

# define transformation for the unit interval
unit <- function (value = 0.5) {
  
  createTransformConstructor(type = 'unit',
                             value = value,
                             true2cont = function(x) qlogis(x),
                             cont2true = function(x) plogis(x),
                             test = function(x) x < 1 & x > 0)
  
}

# define transformation for positive numbers
pos <- function (value = 1) {
  
  createTransformConstructor(type = 'pos',
                             value = value,
                             true2cont = function(x) log(x),
                             cont2true = function(x) exp(x),
                             test = function(x) x > 0)
  
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

p <- pos()
p()
p(TRUE)
test(p)(0.1)
test(p)(-0.1)

p <- update(p, 0.1)
p()
p <- update(p, -5, cont = TRUE)
p()
p(TRUE)


# p <- unit()
# p()
# p(TRUE)
# test(p)(0.1)
# test(p)(-0.1)
# 
# p <- update(p, 0.1)
# p()
# p <- update(p, -5, cont = TRUE)
# p()
# p(TRUE)
