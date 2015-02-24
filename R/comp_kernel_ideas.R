# play at packing and unpacking compositional kernels
rm(list = ls())
library(gpe)

k1 <- rbf('a')
k2 <- lin(c('a', 'b'))
k3 <- k1 + k1 * k2
k3

getObject <- gpe:::getObject
class(getObject(k3)$kernel1)
str(k3, 1)



# option 1:
# make a setObject function (convert object list into a function)
# write a repackKernel function by parsing the symbolic string

repackKernel <- function (unpacked_kernel) {
  
  # convert the data_list basis function kernels
  # into kernel functions again
  unpacked_kernel$data_list <- lapply(unpacked_kernel$data_list,
                                      setObject)
  
  # now just need to substitute these back into the symbolic string
  # and evaluate it...
  
}
# playing wiht expressions
tmp <- list(A = 10, B = pi)
string <- '(A + B)'

# could dump list as named arguments in the local namespace!
# could at least do that with assign and a loop

?expression




# option 2:
# use relist functions; write a function to listify a compositional
# kernel, then unlist, modify, and relist it



tmp$symbolic_string

A <- exp
B <- log
C <- pmin
BC <- list(B = B, C = C)
ABC <- list(A = A, BC = BC)

flesh <- unlist(ABC)
ABC2 <- relist(flesh, ABC)
ABC2

flesh

# convert a compositional kernel object into a list via recursion
as.list.kernel <- function (kernel, recursive = TRUE) {
  checkKernel(kernel)
  object <- getObject(kernel)
  
  # for now, throw an error the kernel is compositional
  if (object$type %in% c('sum', 'prod', 'kron')) {
    object[2:3] <- lapply(object[2:3], as.list.kernel)
  } else {
    class(object) <- 'kernel_object'
  }
  
  return (object)
  
}

# think this won't work as the type is an element in the list
# it may need to be an attribute instead.
# shouldn't be a problem for the basis kernels, though that
# might also screw up unlist/relist
kernel_list <- as.list.kernel(k3)
kernel_list_flat <- unlist(kernel_list, recursive = FALSE)


str(kernel_list_flat, 2)

args(getObject)

as.kernel.list <- function (list) {
  
}


