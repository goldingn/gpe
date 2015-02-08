# functions for the kernel class -all S3 classes

#' @name kernel
#' @rdname kernel
#'
#' @title kernel object
#' 
#' @description Generic functions associated with the kernel class.
#' @template kcl_kernel
#' 
#' @examples
#'  
#' # construct a kernel with one feature
#' k1 <- rbf('temperature')
#'  
NULL

as.kernel <- function (x) {
  
  # add kernel to x's class
  if(!inherits(x, "kernel")) {
    class(x) <- c("kernel", class(x))
  }
  
  # return the x
  return (x)
  
}

#' @rdname kernel
#' @export
#' @examples
#' 
#' # is it a kernel? 
#' is.kernel(k1)
#'  
is.kernel <- function (x) {
  
  # test whether x is a kernel function
  ans <- inherits(x, "kernel")
  
  # return the answer
  return (ans)
  
}

# print function for kernels - just a basic structure
#' @rdname kernel
#' @param \dots For compatibility with the generic print function, not used.
#' @export
#' @examples
#'  
#' # print the kernel's basic structure
#' print(k1)
#'  
#' # this is also the default action for displaying the object:
#' k1
#'   
print.kernel <- function (x, ...) {
  
  # parse kernel structure into a readable format
  string_simple <- parseKernelStructure(x)
  
  # trim outer parentheses from these
  string_simple <- trimParentheses(string_simple)
  
  # print human-readable kernel structure
  cat(paste0('\t\t',
             string_simple,
             '\n'))
  
}


# summary function for kernels
#' @rdname kernel
#' @param digits The number of digits to display for the kernel parameters
#' @export
#' @examples
#'  
#' # a more detailed summary of the kernel's structure
#' summary(k1)
#'  
summary.kernel <- function (object, digits = getOption("digits")) {
  
  # get the object's name
  name <- deparse(substitute(object))
  
  # unpack the kernel into a linear structure
  kernel_data <- unpackKernel(object)
  
  # get structure in symbolic format
  string_symbolic <- kernel_data$symbolic_string
  
  # trim outer parentheses from these
  string_symbolic <- trimParentheses(string_symbolic)
  
  # print more detailed kernel information
  cat('\nKernel summary\n\n')
  
  # if it's more than a single basis kernel
  if (nchar(string_symbolic) > 1) {
    
    # print the symbolic structure
    cat(paste0('\t\t\t\t',
               name,
               ' = ',
               string_symbolic,
               '\n\n'))
    
    # loop through the basis kernels printing out their details
    for (i in 1:length(kernel_data[[2]])) {
      cat(paste0('Basis kernel ',
                 LETTERS[i],
                 '\n'))
      dumpDetails(kernel_data[[2]][[i]], digits)
      cat('\n')
    }
    
  } else {
    # otherwise for a single basis function
    dumpDetails(kernel_data[[2]][[1]], digits)
  }
  
}
