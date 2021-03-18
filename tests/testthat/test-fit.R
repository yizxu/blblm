#test_that("multiplication works", {
#  x1 = rnorm(10000)
#  x2 = rnorm(10000)
#  y = 2+3*x1+4*x2
#  options(future.rng.onMisuse="ignore")
#  df = data.frame(cbind(y,x1,x2))
#  fit0 = blblm(y~x1*x2, df, m = 3, B = 1000, parallel = F)
#  fit1 = blblm(y~x1*x2, df, m = 3, B = 1000, parallel = T)
#  expect_equal(as.integer(coef(fit0)), as.integer(coef(fit1)), tolerance = 1)
#  expect_equal(as.integer(coef(fit0)), c(2,3,4,0), tolerance = 1)
#})

test_that("multiplication works", {
  x1 = rnorm(10000)
  x2 = rnorm(10000)
  y = 2+3*x1+4*x2
  options(future.rng.onMisuse="ignore")
  df = data.frame(cbind(y,x1,x2))
  fit0 = blblm(y~x1*x2, df, m = 3, B = 1000, parallel = F)
  fit1 = blblm(y~x1*x2, df, m = 3, B = 1000, parallel = T)
  expect_equal(as.numeric(coef(fit0)), as.numeric(coef(fit1)), tolerance = 0.05)
  expect_equal(as.numeric(coef(fit0)), c(2,3,4,0), tolerance = 0.05)
})