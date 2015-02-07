# functions to handle syntactic sugar and overload existing functions

# ~~~
# overload the arithmetic operators for addition, multiplication
# and the kronecker product

`+.kernel` <- function (k1, k2) kernel.comp(k1, k2, 'sum')

`*.kernel` <- function (k1, k2) kernel.comp(k1, k2, 'prod')

`%x%` <- function (X, Y) kronecker(X, Y)

kronecker <- function (X, Y, ...) {
  # replace the present kronecker function with one that can handle kernels
  # unfortunately can't overload %x% directly, or provide a kronecker.kernel
  # function
  
  # get a copy of the base R kronecker function to avoid infinite recursion
  kronecker_original <- get("kronecker", envir=baseenv())
  
  # if both things are kernels, get a kronecker product kernel
  if (is.kernel(X) & is.kernel(Y)) {
    
    ans <- kernel.comp(k1, k2, 'kron')
    
  } else {
    
    # otherwise, use the original function
    ans <- kronecker_original(X, Y, ...)
    
  }
  
  # return the result
  return (ans)
  
}


# similarly need to overload the sum and product operators,
# & handle multiple inputs