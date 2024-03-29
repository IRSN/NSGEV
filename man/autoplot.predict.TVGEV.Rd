% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictTVGEV.R
\name{autoplot.predict.TVGEV}
\alias{autoplot.predict.TVGEV}
\title{Autplot Predict Results for \code{TVGEV}}
\usage{
\method{autoplot}{predict.TVGEV}(object, bw = TRUE, ...)
}
\arguments{
\item{object}{An object with class \code{"predict.TVGEV"} as returned
by the \code{predict} method. This must be a data frame so the
\code{predict} method must be called without using \code{out} or
setting it to the default value \code{"data.frame"}.}

\item{bw}{Logical. Should the plot render in black and white for
printing?}

\item{...}{Not used yet.}
}
\value{
An object with class \code{"gg"}. Can be used with the method
\code{plot}, or equivalently with \code{print}.
}
\description{
Produce a \pkg{ggplot2} plot of the results of \code{predict}
for a \code{TVGEV} object.
}
\note{
This function is intended to work with predictions computed
for a possibly large number of periods (to get smooth curves) but
with only a small number of dates, each appearing in a facet. So
the number of dates is limited to \code{6}. Similarily, the number
of confidence levels can not be \code{> 3}.
}
\examples{
example(TVGEV)
pred <- predict(res2, newdate = c("1960-01-01", "2000-01-01", "2020-01-01"),
                level = c(0.70, 0.95), confintMethod = "delta")
g <- autoplot(pred)
}
