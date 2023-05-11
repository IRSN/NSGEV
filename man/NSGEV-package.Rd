\name{NSGEV-package}
\alias{NSGEV-package}
\docType{package}
\title{

  \packageTitle{NSGEV}

}
\description{

  \packageDescription{NSGEV}

}
\details{

  The DESCRIPTION file:
  \packageDESCRIPTION{NSGEV}
  %% \packageIndices{NSGEV}
  
  This package is dedicated to some Non-Stationary Extreme Value models
  with Generalised Extreme Value (GEV) marginal distributions. Such
  models are often used to investigate climate change from series of
  annual maxima.

  An emphasis is put on \emph{Profile Likelihood} inference on the
  parameters, quantiles or return levels. For that aim constrained
  optimisations are used, relying on the \pkg{nloptr} package.

  The package provides most of the classical S3 methods useful to cope
  with fitted model objects, such as \code{predict}, \code{coef},
  \code{vcov}, \code{confint}, \code{residuals}, \code{logLik}, ...

  \itemize{
    
    \item \code{"NSGEV"} This class was intended to be the main class
    when the package was created. As a major difference with other EV
    packages, the model parameters can be used in several of the GEV
    parameters \code{loc}, \code{scale} and \code{shape}.  A typical
    example is when a scale proportional to the location needs to be
    specified.
    
    \item \code{"TVGEV"} This class is devoted to the special case of
    \emph{Time Varying} (TV) models using only functions of the time as
    covariates.  Specialised methods make sense to plot, predict and
    more.
    
  }

  Most efforts have been put on the \code{"TVGEV"} class. In the future,
  the \code{"NSGEV"} class could be discarded.
 
  
}
\author{

   \packageAuthor{NSGEV}

   Maintainer: \packageMaintainer{NSGEV}
}
\references{
  
  Coles, S. (2001) \emph{Introdution to Statistical Modelling of Extreme
  Values}, Springer.

}
\keyword{package}
%% \seealso{}
% \examples{}
