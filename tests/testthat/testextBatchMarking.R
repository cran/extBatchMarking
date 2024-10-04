
# Load data
testData <- WeatherLoach

# Initial parameter
theta <- c(0, -1)

thet <- c(0.1, 0.1, 7, -1.5)

res1 <- batchMarkHmmLL(par           = theta,
                       data          = testData,
                       covariate_phi = NULL,
                       covariate_p   = NULL,
                       choiceModel   = "model4")


res2  <- batchMarkUnmarkHmmLL(par           = thet,
                              data          = testData,
                              choiceModel   = "model4",
                              covariate_phi = NULL,
                              covariate_p   = NULL,
                              Umax          = 1800,
                              nBins         = 600)


# Test functions
test_that( "Equality", {

  expect_equal(round(res1,2),       round(132.3349, 2))
  expect_equal(round(res2,2),       round(870.0261, 2))

} )
