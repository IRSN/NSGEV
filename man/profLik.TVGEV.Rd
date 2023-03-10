% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/profLik.R
\name{profLik.TVGEV}
\alias{profLik.TVGEV}
\title{Profile-Likelihood Inference for \code{TVGEV} Objects}
\usage{
\method{profLik}{TVGEV}(object, fun, level = 0.95, deriv = TRUE, trace = 0, ...)
}
\arguments{
\item{object}{A \code{TVGEV} object.}

\item{fun}{A function of the parameter vector for which the
profile-likelihood will be carried over. This function must have
the arguments: \code{psi} for the vector of parameters and
\code{object} for the model object; so the function can use the
some of slots of \code{object}. If needed, a wrapper function can
be used use more arguments, see \bold{Details}.}

\item{level}{Level of confidence. Can be of length \code{> 1}.}

\item{deriv}{Logical. If \code{TRUE}, the function \code{fun} is
assumed to provide a gradient vector as an attribute named
\code{"gradient"} of the result. For now \code{deriv} can only be
\code{TRUE}, which implies that \code{fun} \emph{must} compute the
gradient.}

\item{trace}{Level of verbosity; \code{trace = 0} prints nothing.}

\item{...}{Not used yet.}
}
\value{
An array with the value of the function and the
corresponding Lower and Upper end-points for the given confidence
levels. This array has two attributes with names \code{"diagno"}
and \code{"psi"} which both are arrays. The attributes provide
information about the numerical optimisation and the values of the
vector of parameter that maximised or minimised the function
\code{fun}.
}
\description{
Profile-likelihood inference for \code{TVGEV} objects.
}
\details{
Compute the lower and upper end-points of a profile-likelihood
based confidence interval. The (apparently new) method used here
relies on maximising and minimising the function of interest, say
\eqn{\eta(\boldsymbol(\psi)}{\eta(\psi)}, under the constraint
that the log-likelihood is greater than the maximal log-likelihood
minus a positive quantity \eqn{\delta} depending on the confidence
level. This differs from the usual method which relies on an
univariate zero-finding for the profile-likelihood function (minus
a constant). Remind that each evaluation of the profile requires a
\eqn{p-1} dimensional optimisation. As a major advantage, the new
method does not require a re-parameterisation of the model.
}
\note{
For each confidence limit the numerical optimisation may
fail, in which case the limit will be \code{NA}. Using \code{trace
= 1} can be useful to further check the optimisation. The
\code{Optimisation status} value should be \code{3} or \code{4}
for small changes on the parameter or on the objective. On the
other hand, a value of \code{5} indicates that the maximal number
of iterations was reached, which is considered here as a
failure. The \code{Constaint check} value should be small because
the constraint must be active at the optimum. The \code{gradDist}
is the distance between the two directions of the gradient vectors
(objective and constraint). It should be small as well because the
gradients must be colinear at the optimum (Lagrange conditions).
}
\examples{
df <- within(TXMax_Dijon, Date <- as.Date(sprintf("\%4d-01-01", Year)))

## fit a TVGEV model with constant parameters.
res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
              estim = "nloptr")

## define a function of the parameter vector: here the first component
## This is for illustration only since the the result can be obtained
## using the 'confint' method with \code{method = "proflik"}, which
## gives the confidence intervals for each of the parameters.

myfun <- function(psi, object) {
    res <- psi[1]
    grad <- rep(0.0, object$p)
    grad[1] <- 1
    attr(res, "gradient") <- grad
    res
}

pl <- profLik(object = res1, fun = myfun, deriv = TRUE)

confint(res1, method = "proflik")
}
\references{
Deville Y. (2017) "Profile-likelihood using constrained
optimisation". Unpublished Tech. Report.
}
\author{
Yves Deville
}
