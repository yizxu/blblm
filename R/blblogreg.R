### blblm and blblogreg shared some of the function (ex. map_mean
#' @name blblogreg
#' @title Logistic Regression with Little Bag of Bootstrap
#' @description `blblogreg` is used to do a logistic regression with
#' implementing little bag of bootstrap. The default setting is in `usage`
#' @usage
#' blblogreg(formula, family, data, m = 10, B = 5000, parallel = FALSE)
#' @param formula The formula for linear regression (response ~ terms)
#' @param family A family (function) used for `glm`
#' @param data The data use for linear regression
#' @param m Number of subsample (split data into m groups)
#' @param B Number of bootstrap
#' @param parallel Determine if using parallelization (logical)
#' @return blblogreg: regression result
#' @export
#' @examples
#' blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, TRUE)

blblogreg <- function(formula, family, data, m = 10, B = 5000, parallel = FALSE) {
  data_list <- split_data(data, m)
  if (parallel == TRUE){
    estimates <- future_map(
      data_list,
      ~ glm_each_subsample(formula = formula, family = family, data = ., n = nrow(data), B = B))
  }
  else if (parallel == FALSE){
    estimates <- map(
      data_list,
      ~ glm_each_subsample(formula = formula, family = family, data = ., n = nrow(data), B = B))
  }
  else {
    warning("The parallel parameter has to be logical (TRUE/FALSE)")
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblogreg"
  invisible(res)

}


#' compute the estimates
#' @param formula Model formula for regression
#' @param family A family (function) used for `glm`
#' @param data Data used to do regression
#' @param n Number of rows in `data`
#' @param B Number of bootstrap

glm_each_subsample <- function(formula, family, data, n, B) {
  replicate(B, glm1(formula, family, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @param formula Model formula for regression
#' @param family A family (function) used for `glm`
#' @param data Data used to do regression
#' @param n Number of rows in `data`

glm1 <- function(formula, family, data, n) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  fit <- glm(formula, data, family = family, weights = freqs)
  list(coef = blbcoef(fit), sigma = sigma(fit))
}


#' @name print.blblogreg
#' @title Print the model of blblogreg
#' @param x The result from `blblogreg`
#' @param ... Additional arguments to be passed to the function
#' @return The formula for the linear regression in `blblogreg`
#' @export
#' @method print blblogreg
#' @examples
#'fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' print(fit)

print.blblogreg <- function(x, ...) {
  cat("blblogreg model:", capture.output(x$formula))
  cat("\n")
}

#' @name sigma.blblogreg
#' @title Sigma value of regression
#' @param object The result from `blblogreg`
#' @param confidence determine if the confience interval is needed. (logical, default is F)
#' @param level The significant level (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @return The sigma value(or values from confidence interval)
#' @export
#' @method sigma blblogreg
#' @examples
#' fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' sigma(fit)
#' sigma(fit, confidence = TRUE)

sigma.blblogreg <- function(object, confidence = FALSE, level = 0.95, ...) {
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

#' @name coef.blblogreg
#' @title Coefficients from blblogreg
#' @param object The result from `blblogreg`
#' @param ... Additional arguments to be passed to the function
#' @return The regression coefficient result from `blblogreg`
#' @export
#' @method coef blblogreg
#' @examples
#' fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' coef(fit)

coef.blblogreg <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @name confint.blblogreg
#' @title Confidence Interval from blblogreg
#' @param object The result from `blblogreg`
#' @param parm The Object (or objects) from `blblogreg` model
#' @param level The significant value. (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @export
#' @method confint blblogreg
#' @examples
#' fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' confint(fit)
#' confint(fit, c("Sepal.Length", "Petal.Length"))

confint.blblogreg <- function(object, parm = NULL, level = 0.95, ...) {
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


#' @name predict.blblogreg
#' @title Prediction for new data with blblogreg
#' @param object The result from `blblogreg`
#' @param new_data New observation that want to predict
#' @param confidence Determine if the confience interval is needed. (logical, default is F)
#' @param level The significant value. (default is 0.95)
#' @param ... Additional arguments to be passed to the function
#' @return The prediction of new data using `blblogreg`
#' @export
#' @method predict blblogreg
#' @examples
#' fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
#' predict(fit, data.frame(Sepal.Length = c(5.0, 5.1), Petal.Length = c(1.3, 1.4)))
#' predict(fit, data.frame(Sepal.Length = c(5.0, 5.1), Petal.Length = c(1.3, 1.4)), confidence = TRUE)

predict.blblogreg <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
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



