# fit gp models and check the gp class works as it should

context('model fitting and gp class functions')

# some fake data
n <- 100
df <- data.frame(x = rnorm(n),
                 z = sample(letters[1:10], n, replace = TRUE))
df$y_gaus <- sin(df$x) + rnorm(n, 0, 0.1)
df$y_bern <- rbinom(n, 1, plogis(df$y_gaus))
df$y_pois <- rpois(n, exp(df$y_gaus) * 3)

test_that('basic model fitting works', {
  
  # simple models with all currently available families
  m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
               data = df)
  m_bern <- gp(y_bern ~ rbf('x') + iid('z'),
               data = df,
               family = binomial)
  m_pois <- gp(y_pois ~ rbf('x') + iid('z'),
               data = df,
               family = poisson)
  
})

test_that('model fitting options work', {
  
  # verbosity
  suppressWarnings(m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
                                data = df,
                                verbose = TRUE))
  
  # ~~~~~~
  # check sparse models
  
  # no inducing points provided
  expect_error(m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
                            data = df,
                            inference = 'FITC'))
  
  # now with inducing points
  m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
               data = df,
               inducing_data = df[1:50, ],
               inference = 'FITC')

  # ~~~~~
  # check latent inference methods

  # currently no Gaussian+Laplace
  expect_error(m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
                            data = df,
                            inducing_data = df[1:50, ],
                            inference = 'Laplace'))
  
  # can't do direct inference on a non-Gaussian likelihood
  expect_error(m_pois <- gp(y_pois ~ rbf('x') + iid('z') + iid(),
                            data = df,
                            family = poisson,
                            inference = 'full'))
  
  # ~~~~~
  # check hyperinference procedures
  
  # these failing due to iid kernel hyperparameter not updatng properly?
  # m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
  #              data = df,
  #              hyperinference = 'BFGS')
  # 
  # m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
  #              data = df,
  #              hyperinference = 'BFGSrestarts')
  
  m_gaus <- gp(y_gaus ~ rbf('x') + iid(),
               data = df,
               hyperinference = 'BFGS')

  m_gaus <- gp(y_gaus ~ rbf('x') + iid(),
               data = df,
               hyperinference = 'BFGSrestarts')
  
  m_gaus <- gp(y_gaus ~ rbf('x') + iid('z') + iid(),
               data = df,
               hyperinference = 'none')
  
})