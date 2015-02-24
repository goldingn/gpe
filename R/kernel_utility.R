# utility functions

trimParentheses <- function (string) {
  # check whether the string starts with a parenthesis and remove outer
  # parentheses if so.
  
  if (substr(string, 1, 1) == '(') {
    string <- substr(string,
                     2,
                     nchar(string) - 1)
  }
  
  return (string)
  
}

dumpDetails <- function(x, digits) {
  # prettily-print a basis kernel's details.
  # x is the list returned by getObject for a basis kernel
  
  cat('\t\t\t\ttype: ',
      x$type)
  cat('\n\t\t\t\tactive columns: ',
      paste(x$columns, collapse = ', '))
  cat('\n\t\t\t\tparameters: ')
  for(i in 1:length(x$parameters)) {
    cat(names(x$parameters)[i],
        ' = ',
        paste(round(x$parameters[[i]], digits), collapse = ', '),
        '\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t')
  }
}

parseKernelStructure <- function (x,
                                  type = c('simple',
                                           'call')) {
  # function to parse structure from nested kernels ready for printing
  # if type = 'simple' then the text returned is a simplified, human-readable
  # structure, e.g. "(rbf(a, b) * rbf(c))"; if (type = 'call' then it
  # is roughly what is needed to recreate the kernel structure, e.g.
  # "(rbf(c('a', 'b')) * rbf('c'))".
  
  
  # match the parsing type
  type <- match.arg(type)
  
  # get the kernel's data object
  x_obj <- getObject(x)
  
  if (x_obj$type %in% c('prod', 'sum', 'kron')) {
    
    # if it's a compositional kernel, recursively call
    # this function to parse subkernels
    string1 <- parseKernelStructure(x_obj$kernel1,
                                    type = type)
    
    string2 <- parseKernelStructure(x_obj$kernel2,
                                    type = type)
    
    # combine the strings, and add parentheses
    string <- paste0('(',
                     string1,
                     switch(x_obj$type,
                            sum = ' + ',
                            prod = ' * ',
                            kron = ' %x% '),
                     string2,
                     ')')
  } else {
    # otherwise, if it's a regular kernel
    
    # if it's wanted in call form
    if (type == 'call') {
      
      # combine the column names, keeping the quotations
      
      # add inverted commas
      columns <- paste0("'",
                        x_obj$columns,
                        "'")
      
      # if more than one column
      if (length(x_obj) > 1) {
        
        # combine with commas
        columns <- paste0(columns, collapse = ', ')
        
        # and add concatenation text
        columns <- paste0('c(', columns, ')')
        
      }
      
    } else if (type == 'simple') {
      
      # otherwise in slightly simplified form
      columns <- paste0(x_obj$columns, collapse = ', ')
      
    }
    
    
    # append the kernel type
    string <- paste0(x_obj$type, '(', columns, ')')
    
  }
  
  # return the resulting string
  return (string)
  
}


# unpacking compositional kernel structure for reporting parameters etc.
# to the user

# unpacked structure is a list with contents:
# symbolic:
#   a text string, giving a symbolic representation of kernel structure, e.g.
#   '((A + B) * C)'
# data:
#   a named list giving the kernel types, columns and parameters of the basis
#   kernels corresponding to the symbols in `symbolic`, e.g.
#     A = list(type = 'rbf',
#              columns = c('a', 'b'),
#              parameters = list(sigma = 1,
#                                l = c(1, 1))),
#     B = ...

unpackKernel <- function (kernel) {
  
  # call flatten kernel
  kernel_data <- flattenKernel(kernel)
  
  # rename the elements of data_list to the letters in symbolic_string
  names(kernel_data$data_list) <- LETTERS[1:kernel_data$counter]
  
  # return the string and data
  return (kernel_data[1:2])
  
}

# do the opposite, take an unpacked kernel and repack it
repackKernel <- function (unpacked_kernel) {
  
  # convert the data_list basis function kernels
  # into kernel functions again
  unpacked_kernel$data_list <- lapply(unpacked_kernel$data_list,
                                      setObject)
  
  # now substitute these back into the symbolic string and evaluate
  kernel <- eval(parse(text = unpacked_kernel$symbolic_string),
                 envir = unpacked_kernel$data_list)
  
  return (kernel)
  
}

# build the following two into getParameters / setParameters
# to enable users to change parameters after constructing kernels

# get the parameters of a kernel (or unpacked kernel)
# as a vector, or vector-like list
getParVec <- function (x) {
  
  # if x is a kernel, unpack it first
  if (is.kernel(x)) {
    x <- unpackKernel(x)
  }
  
  # get parameters only as an unstructured list or vector
  par_list <- lapply(x$data_list,
                     function (x) x$parameters)
  
  # and as a vector or flat list
  par_vec <- unlist(par_list)
  
  return (par_vec)
  
}

# given a kernel (or unpacked kernel) and a corresponding vector
# or vector-like list, such as one given by getParVec, replace the
# the parameters with the elements of this list
setParVec <- function (x, par_vec) {
  
  # flag for whether to return as a kernel or an unpacked kernel
  was_kernel <- FALSE
  
  # if x is a kernel, unpack it first
  if (is.kernel(x)) {
    was_kernel <- TRUE
    x <- unpackKernel(x)
  }
  
  # get parameters only as an unstructured list or vector
  par_list <- lapply(x$data_list,
                     function (x) x$parameters)
  
  # relist the vector into this structure
  par_list <- relist(par_vec, par_list)
  
  # loop through basis kernels, relpacing the parameters
  for (basis in 1:length(par_list)) {
    x$data_list[[basis]]$parameters <- par_list[[basis]]
  }
  
  if (was_kernel) {
    x <- repackKernel(x)
  }
  
  return (x)
  
}

# function to recurse through a kernel object and flatten the structure
flattenKernel <- function (kernel, counter = NULL, data_list = NULL) {
  # counter is used to keep track of the components so far
  # and data_list stores the basis finction details
  
  # if at the top, initialise the counter
  if (is.null(counter)) counter <- 0
  
  # get the kernel's data object
  x_obj <- getObject(kernel)
  
  if (x_obj$type %in% c('prod', 'sum', 'kron')) {
    # if a compositional kernel, unpack the two sub-kernel, passing the
    # counter and data_list back and forth
    
    part_one_list <- flattenKernel(x_obj$kernel1, counter, data_list)
    string1 <- part_one_list[[1]]
    data_list <- part_one_list[[2]]
    counter <- part_one_list[[3]]
    
    part_two_list <- flattenKernel(x_obj$kernel2, counter, data_list)
    string2 <- part_two_list[[1]]
    data_list <- part_two_list[[2]]
    counter <- part_two_list[[3]]
    
    # combine the symbolic strings, and add parentheses
    symbolic_string <- paste0('(',
                              string1,
                              switch(x_obj$type,
                                     sum = ' + ',
                                     prod = ' * ',
                                     kron = ' %x% '),
                              string2,
                              ')')
    
  } else {
    
    # otherwise, if it's a basis kernel
    
    # increment the counter
    counter <- counter + 1
    
    # stop it form re-using letters
    if (counter > 26) {
      
      stop ('cannot currently represent kernels with more than 26 basis kernels symbolicly')
      
    }
    
    # get the next available letter as the symbolic string
    symbolic_string <- LETTERS[[counter]]
    
    # add the kernel object data to data_list
    data_list[[counter]] <- x_obj
    
  }
  
  # return as a list
  ans <- list(symbolic_string = symbolic_string,
              data_list = data_list,
              counter = counter)
  
  return (ans)
  
}


rgp <- function(n, kernel, data = NULL) {
  # draw n random GPs from a kernel, evaluated against a data.frame
  
  # if data is not provided  
  if (is.null(data)) {
    
    # make some new data
    data <- getFakeData(kernel)    
    
  } else {
    
    # otherwise, check it's a dataframe
    stopifnot(class(data) == 'data.frame')
    
  }
  
  # evaluate the kernel on this data
  K <- kernel(data)
  
  # get the number of datapoints
  N <- nrow(data)
  
  # get a random matrix of the correct size to multiply against
  z <- matrix(rnorm(n * N),
              nrow = N,
              ncol = n)
  
  # get random draws
  ans <- K %*% z
  
  # name the columns
  colnames(ans) <- paste0('draw_', 1:n)
  
  # return the draws
  return (ans)
  
}


getFakeData <- function (kernel,
                         n = 1000) {
  # create fake data matching the column names of a kernel
  
  # for now, throw an error the kernel is compositional
  if (getType(kernel) %in% c('sum', 'prod', 'kron')) {
    stop ("sorry, this functionality is currently only available for basis kernels")
  }
  
  # get the column names and dimension
  columns <- getColumns(kernel)
  D <- length(columns)
  
  # random uniform draws for each column
  ans <- replicate(D, runif(1000, -5, 5), simplify = FALSE)
  
  # add their names and turn into a dataframe
  names(ans) <- columns
  ans <- as.data.frame(ans)  
  
  return (ans)

}

expandFactor <- function (factor) {
  
  # stick the factor in a dataframe
  df <- data.frame(x = factor)
  
  # convert into full contrasts matrix
  ans <- model.matrix(~ x - 1, data = df)
  
  return (ans)
  
}

# basis kernel meta-constructor

# This returns a kernel function (R function of class kernel)
# with the kernel type `type`, active on columns `columns` and with parameter
# list `parameters`. The aim of this function is to facilitate kernel function
# creation and reduce the amount of duplicated code
createKernelConstructor <- function(type,
                                    columns,
                                    parameters,
                                    evaluator) {
  
  # create the model object and initialize parameters
  object <- list(type = type,
                 columns = columns,
                 parameters = parameters)
  
  # create a function to return
  ans <- function(data,
                  newdata = NULL,
                  diag = FALSE) {
    
    evaluator(object,
              data,
              newdata,
              diag = diag)    
    
  }
  
  # tell this function it's now a kernel
  ans <- as.kernel(ans)
  
  return(ans)
  
}


# check that the kernel isn't being asked to evaluate the diagonal
# of a non-symmetric matrix
checkSymmetric <- function(newdata) {
  if (!is.null(newdata)) {
    stop ("can't evaluate a diagonal matrix unless newdata is NULL")
  }
}

# diagonal kernel
diagSigma <- function(object, data) {
  # if the diagonal of kernel(data) is just a parameter sigma
  # squared times identity, and object is the kernel's object,
  # return this diagonal matrix - saves code when writing the diagonal
  # option of a kernel evaluator
  ans <- diag(nrow(data)) * object$parameters$sigma ^ 2
}


