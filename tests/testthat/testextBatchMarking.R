
# Load data
testData <- WeatherLoach

# Initial parameter
theta <- c(0, -1)

cores <- detectCores()-1

thet <- c(0.1, 0.1, 7, -1.5)

res1 <- batchMarkHmmLL(par         = theta,
                       data        = testData,
                       choiceModel = "model4",
                       cores)


res2  <- batchMarkUnmarkHmmLL(par         = thet,
                              data        = testData,
                              choiceModel = "model4",
                              Umax        = 1800,
                              nBins       = 20,
                              cores       = cores)


# Test functions
test_that( "Equality", {

  expect_equal(round(res1,2),       round(132.3349, 2))
  expect_equal(round(res2,2),       round(394.3386, 2))

} )
