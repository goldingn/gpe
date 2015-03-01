# functions for the kernel class

#' @name kernel
#' @rdname kernel
#'
#' @title gpe kernel class
#' 
#' @description Generic functions associated with the kernel class.
#' @template kcl_kernel
#' 
#' @details \code{is.kernel} returns a logical indicating whether the 
#' object is a gpe kernel. \code{print} returns a very simple summary of the 
#' kernel structure. \code{summary} returns a more detailed summary of the 
#' kernel structure, including the values of the kernel parameters.
#' \code{plot} plots the covariance structure of the kernel (how covariance)
#' between two points depends on the distance between them.
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
  ans <- inherits(x, "kernel") &
    inherits(x, "function")
  
  # return the answer
  return (ans)
  
}

# print function for kernels - just a basic structure
#' @rdname kernel
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
#' @param object a kernel object
#' @param digits the number of digits to display for the kernel parameters
#' @export
#' @examples
#'  
#' # a more detailed summary of the kernel's structure
#' summary(k1)
#'  
summary.kernel <- function (object, ..., digits = max(3, getOption("digits")-3)) {
  
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



#' @rdname kernel
#' @param which a string giving the column name for which to plot the 
#' covariance structure. If \code{NULL} a plot will be produced for each 
#' caolumn on which the kernel acts
#' @param range the range of distances between points over which to plot 
#' covariances
#' @param ref the reference value against which to compare other points 
#' of a given distance
#' 
#' @export
#' @examples
#'  
#' # plot the covariance structure of the GP
#' plot(k1)
#' 
plot.kernel <- function(x, which = NULL, range = c(-5, 5), ref = 0, ...) {
  
  # check x is a kernel object
  checkKernel(x)
  
  # check ref and range are sensible
  if (! ((ref <= range[2] & ref >= range[1]) |
           (ref >= range[2] & ref <= range[1]))) {
    stop ('ref must fall between the two values of range')
  }
  
  stopifnot(is.finite(ref))
  stopifnot(all(is.finite(range)))
  # get the columns on which the kernel acts
  columns <- getColumns(x)
  
  # number of columns
  D <- length(columns)
  
  if (is.null(which)) {
    # if they want all columns plotted, recursively call this function
    
    for (column in columns) {
      plot(x,
           which = column,
           range = range,
           ref = ref,
           ...)
      
      title(main = paste0('Covariance of the kernel:\n',
                          trimParentheses(parseKernelStructure(x)),
                          '\nagainst ', column))
    }
    
    
    
  } else {
    
    # otherwise act only on that column
    
    # check it's a valid column
    if (!(which %in% columns)) {
      stop (paste0('the kernel does not act on column ',
                   which,
                   ', it only acts on columns: ',
                   paste(columns, collapse = ', ')))
    }
    
    # number of points in the line to predict to
    n_predict <- 500
    
    # create dummy dataframes containing 0s for all variables
    df_from <- data.frame(replicate(D,
                                    0,
                                    simplify = FALSE))
    
    df_to <- data.frame(replicate(D,
                                  rep(0, n_predict),
                                  simplify = FALSE))
    
    # give them the correct names
    names(df_from) <- names(df_to) <- columns
    
    # assign to the column of interest
    df_from[, which] <- ref
    df_to[, which] <- seq(range[1], range[2], len = 500)
    
    # evaluate the kernel on these dataframes
    K <- x(df_from, df_to)
    
    plot(K[1, ] ~ df_to[, which],
         type = 'l',
         lwd = 3,
         col = grey(0.3),
         ylab = 'K_ij',
         xlab = paste0(which,
                       '_i - ',
                       which,
                       '_j'),
         ...)
  }
}