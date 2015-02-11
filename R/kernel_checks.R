# kernel_checks


checkParameterNames <- function (kernel, parameter_list) {
  # check the names of new_parameters check those in current_parameters
  # and issue a helpful warning if not
  
  # make sure there are no duplicated names in the new parameter list
  if (any(duplicated(names(parameter_list)))) {
    
    stop('there appear to be duplicated entries in the new parameter list')
    
  }
  
  # extract current parameters
  current_parameters <- getParameters(kernel)
  
  # check the names match
  name_match <- names(parameter_list) %in% names(current_parameters)
  
  # if there's one which doesn't, throw an error
  if (!all(name_match)) {
    
    # find the mismatch
    mismatch <- names(parameter_list)[which(!name_match)]
    if (length(mismatch) == 1) {
      stop (paste0(mismatch,
                   ' is not a valid parameter of a ',
                   getType(kernel),
                   ' kernel.\nValid parameters are: ',
                   paste0(names(current_parameters),
                          collapse = ', ')))
    } else if (length(mismatch) > 1) {
      stop (paste0(paste0(mismatch,
                          collapse = ', '),
                   ' are not valid parameters of a ',
                   getType(kernel),
                   ' kernel.\nValid parameters are: ',
                   paste0(names(current_parameters),
                          collapse = ', ')))
    }
  }
}

checkCompositional <- function (kernel, not = FALSE) {
  # throw an error if the kernel is (or not) compositional
  
  # get kernel type
  type <- getType(kernel)
  
  # see if it's compositional
  is_compositional <- type %in% c('sum', 'prod', 'kron')
  
  # if it's not supposed to be compositional
  if (not & is_compositional) {
    
    stop(paste0('this function can only be used with basis kernels, ',
                'this is a compositional kernel of type: ',
                type))
    
  } else if (!not & !is_compositional) {
    
    stop(paste0('this function can only be used with compositional kernels, ',
                'this is a basis kernel of type: ',
                type))
  }
  
}

checkParameterType <- function (kernel, parameter_list) {
  # check the parameters are of the right size
  # (and later add domain check)
  
  # get current parameters
  parameters_current <- getParameters(kernel)
  
  # check each element in parameter_list
  for (i in 1:length(parameter_list)) {
    
    # get the new parameter name
    parameter_name <- names(parameter_list)[i]
    
    # and length
    parameter_len <- length(parameter_list[[i]])
    
    # find matching parameter
    j <- match(parameter_name,
               names(parameters_current))
    
    # get length of this target parameter
    target_parameter_len <- length(parameters_current[[j]])
    
    # check it is of the correct dimension
    if (target_parameter_len != parameter_len) {
      stop (paste0(parameter_name,
                   'is of length ',
                   parameter_len,
                   'but should be of length ',
                   target_parameter_len))
    }
  }
}

checkParameters <- function (kernel, parameter_list) {
  # check the parameter list fits the kernel
  
  # make sure it's not compositional
  checkCompositional(kernel, not = TRUE)
  
  # make sure the names match
  checkParameterNames(kernel, parameter_list)
  
  # check they're the right type, length etc.
  checkParameterType(kernel, parameter_list)
  
}

checkKernel <- function (kernel) {
  # check this thing is a kernel
  
  if (!is.kernel(kernel)) {
    stop (paste0('a kernel object is required, but an object of class ',
                 paste0(class(kernel), collapse = ', '),
                 ' was passed instead'))
  }
  
}