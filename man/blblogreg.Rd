% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblogreg.R
\name{blblogreg}
\alias{blblogreg}
\title{Logistic Regression with Little Bag of Bootstrap}
\usage{
blblogreg(formula, family, data, m = 10, B = 5000, parallel = FALSE)
}
\arguments{
\item{formula}{The formula for linear regression (response ~ terms)}

\item{family}{A family (function) used for \code{glm}}

\item{data}{The data use for linear regression}

\item{m}{Number of subsample (split data into m groups)}

\item{B}{Number of bootstrap}

\item{parallel}{Determine if using parallelization (logical)}
}
\value{
blblogreg: regression result
}
\description{
\code{blblogreg} is used to do a logistic regression with
implementing little bag of bootstrap. The default setting is in \code{usage}
}
\examples{
blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, FALSE)
blblogreg(Species ~ Sepal.Length * Petal.Length, binomial, iris[1:100,], 3, 100, TRUE)
}
