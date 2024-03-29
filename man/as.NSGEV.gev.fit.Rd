% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asNSGEV.R
\name{as.NSGEV.gev.fit}
\alias{as.NSGEV.gev.fit}
\title{Coerce a \code{gev.fit} Object into a \code{NSGEV} Object}
\usage{
\method{as.NSGEV}{gev.fit}(x, data, ...)
}
\arguments{
\item{x}{The object to be coerced.}

\item{data}{The matrix of covariates that was used.}

\item{...}{Not used yet.}
}
\value{
An object with class \code{"NSGEV"}.
}
\description{
Coerce a \code{gev.fit} object from the \pkg{ismev} package into a
\code{NSGEV} object.
}
\note{
A \code{gev.fit} object does not embed the set of covariates
used nor a even a call. Thus the set of covariates must be given.
To be compliant with \code{NSGEV}, \emph{this data set must be a data frame
or a matrix with colnames}.
}
\examples{
require(ismev)
## number of observations and covariates
n <- 20; m <- 4 

## generate a matrix covariates
set.seed(1234)
dat <- matrix(runif(n * m), nrow = n)
colnames(dat) <- rNames(m)

## response
y <- drop(nieve::rGEV(n))

## fit 
fit <- gev.fit(xdat = y,             ## response     
               ydat = dat,           ## matrix of covariates
               mul = 1:2, sigl = 4,  ## indices of covariates
               shl = NULL,           
               mulink = exp)         ## optional inverse-link

## coerce
ns <- as.NSGEV(fit, data = dat)

}
\author{
Yves Deville
}
