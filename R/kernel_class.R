# functions for the kernel class -all S3 classes

# coercion and testing for the new kernel class
as.kernel <- function (x) {
  
  # add kernel to x's class
  if(!inherits(x, "kernel")) {
    class(x) <- c("kernel", class(x))
  }
  
  # return the x
  return (x)
  
}

is.kernel <- function (x) {
  
  # test whether x is a kernel function
  ans <- inherits(x, "kernel")
  
  # return the answer
  return (ans)
  
}

# print function for kernels - just a basic structure
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
