test_that("multiplication works", {
  options(future.rng.onMisuse = "ignore")
  options(warning = FALSE)
  fit0 <-
    blblogreg(formula = Species ~ Sepal.Length * Petal.Length,
              family = binomial, data = iris[1:100,], m = 3, B = 100,
              parallel = FALSE)
  fit1 <-
    blblogreg(formula = Species ~ Sepal.Length * Petal.Length,
              family = binomial, data = iris[1:100,], m = 3, B = 100,
              parallel = TRUE)
  expect_equal(as.numeric(coef(fit0)), as.numeric(coef(fit1)), tolerance = 50)
})