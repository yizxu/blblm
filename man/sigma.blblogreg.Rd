% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblogreg.R
\name{sigma.blblogreg}
\alias{sigma.blblogreg}
\title{Sigma value of regression}
\usage{
\method{sigma}{blblogreg}(object, confidence = FALSE, level = 0.95, ...)
}
\arguments{
\item{object}{The result from \code{blblogreg}}

\item{confidence}{determine if the confience interval is needed. (logical, default is F)}

\item{level}{The significant level (default is 0.95)}

\item{...}{Additional arguments to be passed to the function}
}
\value{
The sigma value(or values from confidence interval)
}
\description{
Sigma value of regression
}
\examples{
fit <- blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
sigma(fit)
sigma(fit, confidence = TRUE)
}
