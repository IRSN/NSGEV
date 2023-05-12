
##'
##' 
##' @title Quantile of a Random Maximum
##' 
##' @param object An fitted model such as created by using
##'     \code{\link{TVGEV}}.
##' 
##' @param ... Arguments for methods.
##' 
##' @return The quantile of a random maximum.
##'
##' @export
##' 
quantMax <- function(object, ...) {
    UseMethod("quantMax")
}


##' Compute the quantiles for the random maximum on a given period or
##' collection of blocks of interest.
##'
##' Let \eqn{M^\star := \max_{b} Y_b}{M = max_b Y_b} be the maximum over
##' the blocks \eqn{b} of interest. Since the blocks are
##' assumed to be independent the distribution function of \eqn{M^\star}
##' is given by
##' \deqn{
##'    F_{M^\star}(m^\star; \, \boldsymbol{\psi}) = \prod_{b} F_{\texttt{GEV}}(m^\star; \,
##'    \boldsymbol{\theta}_b)
##' }{F_M(m; psi) = prod_b F_GEV(m; \theta_b)}
##' and it depends on the vector
##' \eqn{\boldsymbol{\psi}}{\psi} of the model parameters
##' through
##' \eqn{\boldsymbol{\theta}_b(\boldsymbol{\psi})}{\theta_b(\psi)}.
##' For a given probability \eqn{p}, the corresponding quantile
##' \eqn{q_{M^\star}(p;\,\boldsymbol{\psi})}{q_M(p; \psi)} is the
##' solution \eqn{m^\star}{m} of
##' \eqn{F^\star(m^\star;\,\boldsymbol{\psi}) = p}{F_M(m; \psi) = p}.
##' The derivative of the quantile w.r.t. \eqn{\boldsymbol{\psi}}{\psi}
##' can be obtained by the implicit function theorem and then be used for
##' the inference e.g., using the "delta method".
##' 
##'
##' @title Quantiles of the Maximum on a Period of Time
##' 
##' @param object A \code{TVGEV} model object.
##' 
##' @param prob A numeric vector giving the (non-exceedance)
##'     probabilities at which the quantiles will be computed.
##'
##' @param date A vector that can be coerced to the class
##'     \code{"Date"} giving the beginnings of the blocks for a period
##'     of interest.
##' 
##' @param level The confidence level.
##' 
##' @param psi Optional vector of coefficients. \bold{Caution} not
##'     tested yet.
##' 
##' @param confintMethod Character indicating the method to be used
##'     for the confidence intervals on the quantiles.
##' 
##' @param out Character indicating what type of object will be
##'     returned. When \code{out} is \code{"data.frame"} the output
##'     actually has the (S3) class \code{"quantMax.TVGEV"} inheriting
##'     from \code{"data.frame"}. A few methods exist for this class.
##' 
##' @param trace Integer level of verbosity.
##' 
##' @param ... Not used.
##'
##' @return An object with its class depending on the value of
##'    \code{out}.
##' 
##'    \itemize{
##'         \item {\code{out = "data.frame" }}{
##'             An object inheriting
##'             from the \code{"data.frame"} class which can be used
##'             in methods such as \code{autoplot}.
##'         }
##'         \item {\code{out = "array" }}{
##'             A 3-dimensional array with dimensions, \emph{probability},
##'             \emph{type of result} (Quantile, Lower or Upper confidence
##'             limit) and \emph{confidence level}.
##'         }
##'   } 
##'
##' @method quantMax TVGEV
##'
##' @export 
##'
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' ## fit a TVGEV model. Only the location parameter is TV.
##' res <- TVGEV(data = df, response = "TXMax", date = "Date",
##'              design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'              loc = ~ t1 + t1_1970)
##' qM1 <- quantMax(res, level = c(0.95, 0.70))
##' qM2 <- quantMax(res,
##'                 date = as.Date(sprintf("%4d-01-01", 2025:2055)),
##'                 level = c(0.95, 0.70))
##' qM2
##' gg <- autoplot(qM2, fillConf = TRUE)
##' gg <- gg + ggtitle("Quantile of the maximum over years 2025-2055")
##' gg
##' 
quantMax.TVGEV <- function(object,
                           prob,
                           date = NULL,
                           level = 0.95,
                           psi = NULL,
                           confintMethod = "delta",
                           out = c("data.frame", "array"),
                           trace = 1L,
                           ...) {
    
    DEBUG <- FALSE
    
    out <- match.arg(out)
    confintMethod <- match.arg(confintMethod)
    
    p <- object$p
    parNames.GEV <- c("loc", "scale", "shape")
    
    if (missing(prob)) {
        prob <-  c(seq(from = 0.80, to = 0.99, by = 0.01),
                   0.995, 0.998, 0.999, 1.0 - 1e-4)
    } else if (any(prob <= 0) || any(prob) >= 1.0) {
        stop("'prob' values must be between 0 and 1")
    }
    prob <- sort(prob)
    fProb <- format(prob)
    nProb <- length(prob)
        
    if (is.null(psi)) psi <- object$estimate

    ## take into account the order. 
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
   
    if (!is(object, "TVGEV")) {
        stop("For now, this method only works when `object' has ",
             "class \"TVGEV\"")
    }
    
    if (is.null(date)) date <- object$fDate
    date <- as.Date(date)

    ## =========================================================================
    ## Compute the vectors 'theta' of GEV parameters, stored as the
    ## rows of a matrix with 3 columns.
    ## =========================================================================
    
    theta <- psi2theta(model = object, psi = psi, date = date, deriv = FALSE,
                       checkNames = FALSE)
    
    nd <- nrow(theta)

    if (!all(object$isCst)) {
        L <- modelMatrices.TVGEV(object, date = date)
        X <- L$X
    } else X <- NULL
    
    if (confintMethod == "delta") {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        covPsi <- vcov(object)
        qNorm <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)

        Quant <- array(NA,
                       dim = c(Prob = nProb, Lim = 3L, Level = nLevel),
                       dimnames = list(Prob = fProb,
                                       Type = c("Quant", "L", "U"),
                                       Level = fLevel)) 
        
        ## =====================================================================
        ## Find a grid of quantiles which should cover the probability
        ## range.
        ## =====================================================================
    
        ## qMin <- min(nieve::qGEV(prob[1],
        ##                         loc = theta[ , 1],
        ##                         scale = theta[ , 2],
        ##                         shape = theta[ , 3]))
        ## qMax <- max(nieve::qGEV(prob[length(prob)],
        ##                         loc = theta[ , 1],
        ##                         scale = theta[ , 2],
        ##                         shape = theta[ , 3]))

        qMin <- .qMin.TVGEV(p = prob[1], theta)
        qMax <- .qMax.TVGEV(p = prob[length(prob)], theta)
        
        ## rq <- qMax - qMin
        ## qMin <- qMin - 0.1 * rq
        ## qMax <- qMax + 0.1 * rq
        mStarGrid <- seq(from = qMin, to = qMax, length = 100) 
        FMStarGrid <- rep(0.0, length(mStarGrid))
        
        for (i in seq_along(mStarGrid)) {
            FMStarGrid[i]<- prod(nieve::pGEV(mStarGrid[i],
                                             loc = theta[ , 1],
                                             scale = theta[ , 2],
                                             shape = theta[ , 3]))
        }

        ## =====================================================================
        ## An interpolation is required to get the values at
        ## the 'pretty' probability values.
        ## =====================================================================
        
        qMStar <- approx(x = FMStarGrid, y = mStarGrid, xout = prob)$y
        
        if (DEBUG) {
            plot(prob, qMStar, pch = 16)
            lines(FMStarGrid, mStarGrid)
        }
        
        logFStar <- rep(0.0, length(prob))
        
        gradpsi <- rep(NA, p)
        names(gradpsi) <- object$parNames
        
        for (iProb in seq_along(prob)) {

            ## =================================================================
            ## - 'Fi', 'fi' are numeric vectors with length 'nDate'
            ##   giving the marginal cdf and the pdf at the quantile.
            ##
            ## - 'gradtheta' is a matrix with nDate rows and 3 columns.
            ## =================================================================
            
            Fi <- nieve::pGEV(qMStar[iProb],
                              loc = theta[ , 1],
                              scale = theta[ , 2],
                              shape = theta[ , 3], 
                              deriv = TRUE)

            ## compute the gradient of log(Fi)
            gradtheta <- attr(Fi, "gradient")
            gradtheta <- gradtheta / Fi
            attr(Fi, "gradient") <- NULL
            
            fi <- nieve::dGEV(qMStar[iProb],
                              loc = theta[ , 1],
                              scale = theta[ , 2],
                              shape = theta[ , 3])
            
            logFStar[iProb] <- sum(log(Fi))
            FStar <- exp(logFStar[iProb])
            fStar <- FStar * sum(fi / Fi)
            
            ## First, 'gradpsi' will be the gradient of log(FStar)
            for (nm in parNames.GEV) {
                if (!object$isCst[nm]) {
                    gradpsi[object$ind[[nm]]] <- 
                        t(gradtheta[ , nm, drop = FALSE]) %*% X[[nm]]
                } else {
                    gradpsi[object$ind[[nm]]] <- sum(gradtheta[ , nm])
                }
            }

            ## Now 'gradpsi' becomes the gradient of the quantile by
            ## the implicit function threorem
            
            gradpsi <- - gradpsi * FStar / fStar

            sdQuant <- as.vector(sqrt(t(gradpsi) %*% covPsi %*% gradpsi))
            
            Quant[iProb, "Quant",  ] <- qMStar[iProb] 
            Quant[iProb, "L", ] <- qMStar[iProb] + sdQuant * qNorm[ , 1L] 
            Quant[iProb, "U", ] <- qMStar[iProb] + sdQuant * qNorm[ , 2L]
            
        }

    }
    
    ## =========================================================================
    ## Transform to a data frame. Using the 'data.table' package is
    ## both simple and efficient.
    ## =========================================================================
    
    if (out != "array") {
        L <- list()
        for (iLevel in seq_along(level)) {
            nm <- fLevel[iLevel]
            
            L[[nm]] <- data.frame(Prob = prob,
                             ProbExc = 1.0 - prob,
                             Quant = Quant[ , "Quant", iLevel],
                             L = Quant[ , "L", iLevel],
                             U = Quant[ , "U", iLevel],
                             Level = level[iLevel])
        }
        
        Quant <- data.table::rbindlist(L)
        class(Quant) <- c("quantMax.TVGEV", "data.frame")
    }
    
    Quant
        
}



# *****************************************************************************

##' Find a lower bound for the quantile of the random maximum over a
##' collection of blocks. The distribution function of the maximum is
##' the product of the marginal (GEV) distributions functions for the
##' blocks over which the maximum is defined. By selecting the
##' smallest GEV shape, the smallest scale and the smallest location,
##' we get a GEV distribution which is smaller in stochastic order
##' than all the marginal GEV distributions. We can take the bound as
##' the quantile of this GEV distribution corresponding to the
##' probability \eqn{p^{1/n^\star}}{p^(1/ nStar)}, where
##' \eqn{n^\star}{nStar} is the number of blocks.
##'
##' This function is used by \code{\link{quantMaxFun.TVGEV}} and
##' \code{\link{quantMax.TVGEV}}.
##' 
##' @title Upper Bound for the Quantile of the Random Maximum using a
##'     \code{"TVGEV"} object.
##' 
##' @param theta A numeric matrix containing the GEV "theta"
##'     parameters with one line for each block \eqn{b} . It must have 3
##'     columns corresponding to the GEV parameters \emph{location}
##'     \eqn{\mu_b}, \emph{scale} \eqn{\sigma_b} and \emph{shape}
##'     \eqn{\xi_b} is that order.
##' 
##' @param p A probability.
##' 
##' @keywords internal
##' @export
##' @seealso
##' \code{\link{.qMax.TVGEV}}
.qMin.TVGEV <- function(theta, p) {
    if (ncol(theta) != 3) {
        stop("'theta' must be a numeric matrix with three columns")
    }
    n <- nrow(theta)
    shapeMin <- min(theta[ , "shape"])
    scaleMin <- min(theta[ , "scale"])
    locMin <- min(theta[ , "loc"])
    nieve::qGEV(p = p^(1 / n), loc = locMin, scale = scaleMin, shape = shapeMin)
}

## *****************************************************************************
##' Find an upper bound for the quantile of the random maximum over a
##' collection of blocks. The distribution function of the maximum is
##' the product of the marginal (GEV) distributions functions for the
##' blocks over which the maximum is defined. By selecting the largest
##' GEV shape, the largest scale and the largest location, we get a
##' GEV distribution which is larger in stochastic order than all the
##' marginal GEV distributions. We can take the bound as the quantile
##' of this GEV distribution corresponding to the probability
##' \eqn{p^{1/n^\star}}{p^(1/ nStar)}, where \eqn{n^\star}{nStar} is
##' the number of blocks.
##'
##' ##' This function is used by \code{\link{quantMaxFun.TVGEV}} and
##' \code{\link{quantMax.TVGEV}}.
##' 
##' @title Upper Bound for the Quantile of the Random Maximum using a
##'     \code{"TVGEV"} object.
##' 
##' @param theta A numeric matrix containing the GEV "theta"
##'     parameters with one line for each block \eqn{b} . It must have 3
##'     columns corresponding to the GEV parameters \emph{location}
##'     \eqn{\mu_b}, \emph{scale} \eqn{\sigma_b} and \emph{shape}
##'     \eqn{\xi_b} is that order.
##' 
##' @param p A probability.
##'
##' @keywords internal
##' @export
##' @seealso
##' \code{\link{.qMin.TVGEV}}
.qMax.TVGEV <- function(theta, p) {
    if (ncol(theta) != 3) {
        stop("'theta' must be a numeric matrix with three columns")
    }
    n <- nrow(theta)
    shapeMax <- max(theta[ , "shape"])
    scaleMax <- max(theta[ , "scale"])
    locMax <- max(theta[ , "loc"])
    nieve::qGEV(p = p^(1 / n),
                loc = locMax, scale = scaleMax, shape = shapeMax)
}

##' @title Distribution Function for a Random Maximum
##'
##' @param object An object that can be used to compute/predict the
##'     distribution of the random maximum.
##'
##' @param ... Further arguments for methods.
##' 
##' @return A function (or closure) that can be used to compute the
##'     probability of non-exceedance for arbitraty given
##'     quantiles.
##'
##' @export
##' 
cdfMaxFun <- function(object, ...) {
    UseMethod("cdfMaxFun")
}

## *****************************************************************************
##' Given an object with class \code{"TVGEV"} and a collection of time
##' blocks (or period) defined by \code{date}, the distribution of the
##' random maximum over the period \eqn{M^\star := \max_b Y_b}{MStar =
##' max_b Y_b} is known. The corresponding distribution function can
##' be obtained.
##' 
##' @title Distribution Function for the Random Maximum on a given
##'     Period
##' 
##' @param object An object with class \code{"TVGEV"}.
##' 
##' @param date An object that can be coerced to the class
##'     \code{"Date"}. If not provided this will be taken as the date
##'     attached to \code{object}.
##'
##' @param psi An optional vector of coefficients for
##'     \code{object}. By default the ML estimate as returned by
##'     \code{coef(object)} is used.
##'
##' @param theta An optional matrix with three columns containing GEV
##'     parameters. The colums are in the order \emph{location},
##'     \emph{scale} and \emph{shape}. When this argument is used
##'     neither \code{date} not \code{psi} can be used.
##' 
##' @param ... Not used.
##' 
##' @return A function (more precisely, a closure). This function has
##'     a single formal argument \code{q} representing a quantile, and
##'     it returns the corresponding probability
##'     \eqn{\text{Pr}\{M^\star \leq q \}}{Pr[MStar <= q]}.
##'
##' @method cdfMaxFun TVGEV
##' @export
##'
##' @seealso \code{\link{quantMaxFun.TVGEV}} for the corresponding
##'     quantile function (or closure).
##'
##' @section Caution: When \code{theta} is given \code{model} is not
##'     used. The distribution function is simply that of the maximum
##'     of independent r.vs following GEV distributions with their
##'     parameters given by (the rows of) \code{theta}.
##'
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' ## fit a TVGEV model. Only the location parameter is TV.
##' res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'               design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'               loc = ~ t1 + t1_1970)
##' cdf <- cdfMaxFun(res1)
##' cdf(c(39.0, 41.0))
##' 
##' ## a 'new period'
##' date <- as.Date(sprintf("%4d-01-01", 2025:2055))
##' cdfNew <- cdfMaxFun(res1, date = date)
##' cdfNew(c(39.0, 41.0))
##' 
cdfMaxFun.TVGEV <- function(object, date = NULL, psi = NULL,
                            theta = NULL, ...) {
    
    if (is.null(theta)) {
        if (is.null(psi)) psi <- object$estimate
        if (is.null(date)) date <- object$fDate
        date <- as.Date(date)
        theta <- psi2theta(model = object, psi = psi, date = date)
    } else if (!is.null(date) || !is.null(psi)) {
        stop("when 'theta' is given, neither 'date' nor 'psi' can ",
             "be given")
    }
    
    outFun <- function(q) {
        res <- rep(0.0, length(q))
        for (i in seq_along(q)) {
            res[i] <- prod(nieve::pGEV(q[i],
                                       loc = theta[ , 1],
                                       scale = theta[ , 2],
                                       shape = theta[ , 3]))
        }
        res
    }

    outFun

}

##' @title Quantile Function for a Random Maximum
##'
##' @param object An object that can be used to compute/predict the
##'     distribution of the random maximum.
##'
##' @param ... Further arguments for methods.
##' 
##' @return A function (or closure) that can be used to compute the
##'     quantiles for arbitraty given probabilities.
##'
##' @export
##' 
quantMaxFun <- function(object, ...) {
    UseMethod("quantMaxFun")
}

## *****************************************************************************
##' 
##' Given an object with class \code{"TVGEV"} and a collection of time
##' blocks (or \emph{period}) defined by \code{date}, the distribution
##' of the random maximum over the period \eqn{M^\star := \max_b
##' Y_b}{MStar = max_b Y_b} is known. The corresponding quantile
##' function \eqn{q_{M^\star}(m^\star)}{q_MStar(mStar)} can be
##' obtained.
##' 
##' @title Quantile Function for the Random Maximum on a given
##'     Period
##' 
##' @param object An object with class \code{"TVGEV"}.
##' 
##' @param date An object that can be coerced to the class
##'     \code{"Date"}. If not provided this will be taken as the date
##'     attached to \code{object}.
##'
##' @param psi An optional vector of coefficients for
##'     \code{object}. By default the ML estimate as returned by
##'     \code{coef(object)} is used.
##'
##' @param theta An optional matrix with three columns containing GEV
##'     parameters. The colums are in the order \emph{location},
##'     \emph{scale} and \emph{shape}. When this argument is used
##'     neither \code{date} not \code{psi} can be used.
##' 
##' @param ... Not used.
##' 
##' @return A function (more precisely, a closure). This function has
##'     a single formal argument \code{p} representing a probability
##'     \eqn{p}, and it returns the corresponding quantile \eqn{q}
##'     such that
##'     \eqn{\text{Pr}\{M^\star \leq q \} = p}{Pr[MStar <= q] = p} where
##'     \eqn{M^\star}{MStar} is the random maximum.
##'
##' @method quantMaxFun TVGEV
##' 
##' @export
##'
##' @section Caution: When \code{theta} is given \code{model} is not
##'     used. The quantile function is simply that of the maximum of
##'     independent r.vs following GEV distributions with their
##'     parameters given by (the rows of) \code{theta}.
##' 
##' @seealso \code{\link{cdfMaxFun}} for the corresponding
##'     quantile function (or closure).
##'
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' ## fit a TVGEV model. Only the location parameter is TV.
##' res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'               design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'               loc = ~ t1 + t1_1970)
##' qf <- quantMaxFun(res1)
##' qf(c(0.90, 0.99, 0.999))
##' 
##' ## a 'new period'
##' date <- as.Date(sprintf("%4d-01-01", 2025:2055))
##' qfNew <- quantMaxFun(res1, date = date)
##' qfNew(c(0.90, 0.99, 0.999))
##' 
quantMaxFun.TVGEV <- function(object, date = NULL, psi = NULL, theta = NULL,
                              ...) {

    if (is.null(theta)) {
        if (is.null(psi)) psi <- object$estimate
        if (is.null(date)) date <- object$fDate
        date <- as.Date(date)
        theta <- psi2theta(model = object, psi = psi, date = date)
    } else if (!is.null(date) || !is.null(psi)) {
        stop("when 'theta' is given, neither 'date' nor 'psi' can ",
             "be given")
    }
    
    FMax <- cdfMaxFun.TVGEV(object = object, theta = theta)

    FMaxZero <- function(q, prob) FMax(q) - prob
    
    outFun <- function(probs) {
        res <- rep(0.0, length(probs))
        for (i in seq_along(probs)) {
            interv <- c(.qMin.TVGEV(theta = theta, p = probs[i]),
                        .qMax.TVGEV(theta = theta, p = probs[i]))
            resTry <- try(uniroot(f = FMaxZero,  interval = interv,
                                  prob = probs[i]))
            if (!inherits(resTry, "try-error")) {
                res[i] <- resTry$root
            }
        }
        res
    }
    outFun

}
