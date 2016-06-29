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
                   "iid('a', sigma = 0.1)")

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
    
    expect_true(all.equal(dim(C_self), c(3, 3)))
    expect_true(all.equal(dim(C_proj), c(3, 2)))
    expect_true(all.equal(dim(C_diag), c(3, 3)))
  }
  

})