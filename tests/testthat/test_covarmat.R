# construct kernels and check that all of the access functions
# do what they oughta

context('covarmat methods')

df <- data.frame(a = 1:5,
                 b = 4:0,
                 c = c(1, 1, 1, 1, 1))

test_that('covarmat functions work', {

  # create a covarmat
  k <- rbf('a') + per('b') * expo('a') + iid()
  C <- k(df)
  
  # make sure it is both a covarmat and a matrix
  expect_true(is.covarmat(C))
  expect_true(is.matrix(C))
  
  # suppress plotting, then make sure image works
  pdf(NULL)
  image(C)
  image(C, axes = FALSE, legend = FALSE)
  dev.off()
  
  # check the colour scheme
  expect_equal(palPuRd(3), c("#D4B9DA", "#E7298A", "#67001F"))
})
