#' @import purrr
#' @import furrr
#' @import stats
#' @import future
#' @import utils

#' @importFrom magrittr %>%
#' @importFrom utils "capture.output"
#' @title Using Little Bag of Bootstraps for Linear Regression
#' @aliases blblm-package
#' @description This package do a linear regression with little bag of bootstraps.
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' @name blblm
#' @title Linear Regression with Little Bag of Bootstrap
#' @description `blblm` is used to do a linear regression with
#' implementing little bag of bootstrap. The default setting is in `usage`
#' @usage
#' blblm(formula, data, m = 10, B = 5000, parallel = FALSE)
#' @param formula The formula for linear regression (response ~ terms)
#' @param data The data use for linear regression
#' @param m Number of subsample (split data into m groups)
#' @param B Number of bootstrap
#' @param parallel Determine if using parallelization (logical)
#' @return blblm: regression result
#' @export
#' @examples
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)

blblm <- function(formula, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (parallel == TRUE){
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  else if (parallel == FALSE){
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  }
  else {
    warning("The parallel parameter has to be logical (TRUE/FALSE)")
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#' @param data Data used to split
#' @param m Split data into m subsamples

split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @param formula Model formula for regression
#' @param data Data used to do regression
#' @param n Number of rows in `data`
#' @param B Number of bootstrap

lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @param X Terms
#' @param y Response
#' @param n Number of rows in data

lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit Fit from `blblm`

blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit Fit from `blblm`

blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @name print.blblm
#' @title Print the model of blblm
#' @param x The result from `blblm`
#' @param ... Additional arguments to be passed to the function
#' @return The formula for the linear regression in `blblm`
#' @export
#' @method print blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
#' print(fit)

print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @name sigma.blblm
#' @title Sigma value of regression
#' @param object The result from `blblm`
#' @param confidence determine if the confience interval is needed. (logical, default is F)
#' @param level The significant level (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @return The sigma value(or values from confidence interval)
#' @export
#' @method sigma blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
#' sigma(fit)
#' sigma(fit, confidence = TRUE)

sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' @name coef.blblm
#' @title Coefficients from blblm
#' @param object The result from `blblm`
#' @param ... Additional arguments to be passed to the function
#' @return The regression coefficient result from `blblm`
#' @export
#' @method coef blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
#' coef(fit)

coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @name confint.blblm
#' @title Confidence Interval from blblm
#' @param object The result from `blblm`
#' @param parm The Object (or objects) from `blblm` model
#' @param level The significant value. (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @export
#' @method confint blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
#' confint(fit)
#' confint(fit, c("wt", "hp"))

confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' @name predict.blblm
#' @title Prediction for new data with blblm
#' @param object The result from `blblm`
#' @param new_data New observation that want to predict
#' @param confidence Determine if the confience interval is needed. (logical, default is F)
#' @param level The significant value. (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @return The prediction of new data using `blblm`
#' @export
#' @method predict blblm
#' @examples
#' fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = FALSE)
#' predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)))
#' predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE)

predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
