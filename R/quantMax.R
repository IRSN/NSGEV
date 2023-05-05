
##'
##' 
##' @title Quantile of a Random Maximum
##' 
##' @param object An fitted model such as crated by \code{\link{TVGEV}}.
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
##' is
##' \deqn{
##'    F_{M^\star}(m^\star) = \prod_{b} F_{\texttt{GEV}}(m^\star; \,
##'    \boldsymbol{\theta}_b)
##' }{F_M(m) = prod_b F_GEV(m; \theta_b)}
##' and this distribution function depnds on the vector
##' \eqn{\boldsymbol{\psi}}{\psi} of the model parameters via
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
##' example(TVGEV)
##' qM1 <- quantMax(res2, level = c(0.95, 0.70))
##' qM2 <- quantMax(res2,
##'                 date = as.Date(sprintf("%4d-01-01", 2025:2055)),
##'                 level = c(0.95, 0.70))
quantMax.TVGEV <- function(object,
                           prob,
                           date = NULL,
                           level = 0.95,
                           psi = NULL,
                           confintMethod = "delta",
                           deriv = FALSE,
                           out = c("data.frame", "array"),
                           trace = 1L) {
    
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
    fLevel <- NSGEV:::formatLevel(level)
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
    
    theta <- psi2theta(model = object, psi = psi, date = date, deriv = deriv,
                       checkNames = FALSE)
    
    nd <- nrow(theta)
    
    ## SOG: Save Our Gradient
    if (deriv) dtheta_dpsi <- attr(theta, "gradient")

    if (!all(object$isCst)) {
        L <- NSGEV:::modelMatrices.TVGEV(object, date = date)
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
    
        qMin <- min(nieve::qGEV(prob[1],
                                loc = theta[ , 1],
                                scale = theta[ , 2],
                                shape = theta[ , 3]))
        qMax <- max(nieve::qGEV(prob[length(prob)],
                                loc = theta[ , 1],
                                scale = theta[ , 2],
                                shape = theta[ , 3]))
        rq <- qMax - qMin
        qMin <- qMin - 0.1 * rq
        qMax <- qMax + 0.1 * rq
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
