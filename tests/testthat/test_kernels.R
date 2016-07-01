# construct kernels and check that all of the access functions
# do what they oughta

context('kernel construction and evaluation')

df <- data.frame(a = 1:3,
                 b = 4:2,
                 c = c(1, 1, 1))

test_that('all the kernel constructors work', {
  
  k <- rbf(c('a', 'b'), sigma = 0.1, l = 0.2)
  k <- expo(c('a', 'b'), sigma = 0.1, l = 0.2)
  k <- lin(c('a', 'b'), sigma = 0.1, c = 0.2)
  k <- mat32(c('a', 'b'), sigma = 0.1, l = 0.2)
  k <- mat52(c('a', 'b'), sigma = 0.1, l = 0.2)
  k <- rq(c('a', 'b'), sigma = 0.1, alpha = 0.2)
  k <- per(c('a', 'b'), sigma = 0.1, l = 0.2, p = 1)
  k <- int(sigma = 0.2)
  k <- iid('a', sigma = 0.1)
  k <- iid(sigma = 0.1)
  
})

test_that('some kernel constructors are fussy', {
  
  # int takes no columns
  expect_error(k <- int('a', sigma = 0.2))
  
  # iid canonly act on one column
  expect_error(k <- iid(c('a', 'b'), sigma = 0.1))
  
})

test_that('kernels produce the right covariance matrices', {
  
  expressions <- c("rbf(c('a', 'b'), sigma = 0.1, l = 0.2)",
                   "expo(c('a', 'b'), sigma = 0.1, l = 0.2)",
                   "lin(c('a', 'b'), sigma = 0.1, c = 0.2)",
                   "mat32(c('a', 'b'), sigma = 0.1, l = 0.2)",
                   "mat52(c('a', 'b'), sigma = 0.1, l = 0.2)",
                   "rq(c('a', 'b'), sigma = 0.1, alpha = 0.2)",
                   "per(c('a', 'b'), sigma = 0.1, l = 0.2, p = 1)",
                   "int(sigma = 0.2)",
                   "iid('a', sigma = 0.1)",
                   "iid(sigma = 0.1)")

  # check matrices
  for (expression in expressions) {
    k <- eval(parse(text = expression))
    
    # three types of matrix
    C_self <- k(df)
    C_proj <- k(df, df[1:2, ])
    C_diag <- k(df, diag = TRUE)
    
    expect_true(is.covarmat(C_self))
    expect_true(is.covarmat(C_proj))
    expect_true(is.covarmat(C_diag))
    
    expect_true(is.matrix(C_self))
    expect_true(is.matrix(C_proj))
    expect_true(is.matrix(C_diag))
    
    expect_equal(dim(C_self), c(3, 3))
    expect_equal(dim(C_proj), c(3, 2))
    expect_equal(dim(C_diag), c(3, 3))
  }
  
})


test_that('kernel access functions work', {
  
  # create kernel functions
  k1 <- rbf(c('a', 'b'), sigma = 0.1, l = 0.2)
  k2 <- int(sigma = 0.2)
  k3 <- iid('a', sigma = 0.1)
  k4 <- per(c('a', 'b'), sigma = 0.1, l = 0.2, p = 1)
  k <- k1 + (k2 * k3 + k4)
  
  # check the column names can be extracted
  expect_equal(getColumns(k1), c('a', 'b'))
  expect_equal(getColumns(k2), NULL)
  expect_equal(getColumns(k3), 'a')
  expect_equal(getColumns(k4), c('a', 'b'))
  expect_equal(getColumns(k), c('a', 'b'))

  # check the parameters can be extracted  
  expect_equal(getParameters(k1), list(sigma = 0.1, l = 0.2))
  expect_equal(getParameters(k2), list(sigma = 0.2))
  expect_equal(getParameters(k3), list(sigma = 0.1))
  expect_equal(getParameters(k4), list(p = 1, l = 0.2, sigma = 0.1))

  # check the parameters can be updated
  k1 <- setParameters(k1, l = 0.3)
  expect_equal(getParameters(k1), list(sigma = 0.1, l = 0.3))
  
  # check the types can be extracted
  expect_equal(getType(k1), 'rbf')
  expect_equal(getType(k2), 'int')
  expect_equal(getType(k3), 'iid')
  expect_equal(getType(k4), 'per')
  expect_equal(getType(k), 'sum')
  
  # check the subkernels can be extracted
  expect_equal(getSubKernel(k, 1), k1)
  expect_equal(getSubKernel(k, 2), (k2 * k3) + k4)
  
})

test_that('kernel summaries work', {
  
  expected <- c("",
                "Kernel summary",
                "",
                "\t\t\t\ttype:  rbf",
                "\t\t\t\tactive columns:  a, b",
                "\t\t\t\tparameters: sigma  =  0.1 ",
                "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tl  =  0.2 ",
                "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
  expect_equal(capture.output(summary(k1)), expected)
  
  # try plotting a kernel
  plot(k1)
  plot(k4)
  
})

