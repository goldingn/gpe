# parameter_class

# member functions for the parameter class, representing (hyper) parameters 
# of a gpe model

#' @name parameter
#' @rdname parameter
#'
#' @title gpe parameter class
#' 
#' @description Functions to create and interact with parameter objects
#' 
#' @param value the value of the parameter, on its true scale.
#' 
#' @param x a parameter object, or an object to be tested as one.
#' 
#' @param object a parameter object, or an object to be tested as one.
#' 
#' @param continuous whether to return (or update) with the true value of the 
#' parameter or with the value of its continuous transform.
#' 
#' @param new_value a new value to assign to the parameter, given either on 
#' the true scale, or on the scale of its continuous transform, depending on
#' \code{continuous}.
#' 
#' @param \dots further arguments passed to or from other methods.
#' 
#' @details The constructor functions \code{pos} and \code{unit} return 
#' parameter objects with a set value. These objects caontain member funcitons
#' the constraints of the parameter (it's support) and define functions for 
#' between the true value and a continuous transformation of it. 
#' 
#' \code{unit} creates a parameter constrained to be on the unit interval 
#' (greater than 0 and less than 1).
#' \code{pos} creates a paramter constrained to be positive.
#' 
#' The value (either the true value or a continuous transformation) of the 
#' parameter can be returned by executing the object.
#' The value can be changed using the \code{update} method.
#' \code{is.parameter} returns a logical indicating whether the 
#' object is a gpe parameter object. \code{print} returns a very simple summary
#' of the parameter, it's transformations and its support.
NULL

#' @rdname parameter
#' @export
#' @examples
#' 
#' # create a parameter on the unit interval
#' # with initial value of 0.9
#' p <- unit(0.9)
#' 
#' # return value (on true scale)
#' p()
#' 
#' # return value on continuous scale
#' p(continuous = TRUE)
#' 
unit <- function (value = 0.5) {
  # define parameter on the unit interval
  
  createParameterConstructor(type = 'unit',
                             value = value,
                             true2cont = function(x) qlogis(x),
                             cont2true = function(x) plogis(x),
                             test = function(x) all(x < 1 & x > 0))
  
}

#' @rdname parameter
#' @export
#' @examples
#' 
#' # create a strictly positive parameter
#' # with initial value of 2.5
#' sigma <- pos(2.5)
#' 
#' # return value (on true scale)
#' sigma()
#' 
#' # return value on continuous scale
#' sigma(continuous = TRUE)
#'   
pos <- function (value = 1) {
  # define transformation for positive numbers  
  createParameterConstructor(type = 'pos',
                             value = value,
                             true2cont = function(x) log(x),
                             cont2true = function(x) exp(x),
                             test = function(x) all(x > 0))
  
}

#' @rdname parameter
#' @export
#' @examples
#' 
#' # are these parameters? 
#' is.parameter(p)
#' is.parameter(sigma)
#'  
is.parameter <- function (x) {
  # test whether x is a parameter
  ans <- inherits(x, 'parameter')
  return (ans)
}



#' @rdname parameter
#' @export
#' @examples
#' 
#' # update on the true scale
#' p <- update(p, 0.6)
#' p()
#' 
#' # update on the continuous scale
#' p <- update(p, -1, continuous = TRUE)
#' p()
#' 
update.parameter <- function(object, new_value, continuous = FALSE) {
  # function to update value of a parameter
  # parameter is the parameter,
  # new_value is the new value
  # cont specifies whether new_value is provided on the scale of the 
  # continuous transformation
  
  # check it's a valid parameter
  checkParameter(object)

  # check the length of new_value
  if (length(object()) != length(new_value)) {
    stop ('new_value is the wrong length')
  }
  
  # if the value is provided on the continuous scale
  if (continuous) {
    
    # switch it to the true scale
    new_value <- cont2true(object)(new_value)
    
  }
  
  # test the new value
  if (!test(object)(new_value)) {
    stop ('new_value is not within the constraints of this parameter')
  }
  
  # update the value
  environment(object)$object$value <- new_value
  
  # return the parameter
  return (object)

}

#' @rdname parameter
#' @export
#' @examples
#' 
#' # print summaries of p and sigma and their member functions
#' p
#' sigma
#' 
print.parameter <- function (x, ...) {
  # print function for parameter objects  
  # first check it's a valid parameter object
  checkParameter(x)
  
  # get type
  type <- environment(x)$object$type
  
  # print all of these things
  cat('parameter object of type:',
      type,
      '\n\tcurrent value:\t\t\t\t\t\t\t',
      paste(x(), collapse = ', '),
      '\n\ttest function:\t\t\t\t\t\t\t',
      fun2string(test(x)),
      '\n\tcont2true function:\t\t',
      fun2string(cont2true(x)),
      '\n\ttrue2cont function:\t\t',
      fun2string(true2cont(x)),
      '\n')
  
}

