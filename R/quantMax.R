
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
##' tv <- TVGEV(data = df, response = "TXMax", date = "Date",
##'             design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'             loc = ~ t1 + t1_1970)
##' qM1 <- quantMax(tv, level = c(0.95, 0.70))
##' date2 <- as.Date(sprintf("%4d-01-01", 2025:2054))
##' qM2 <- quantMax(tv, date = date2, level = c(0.95, 0.70))
##' head(qM2)
##' gg <- autoplot(qM2, fillConf = TRUE)
##' gg <- gg + ggtitle("Quantile of the maximum over years 2025-2054 with \"delta\" intervals")
##' gg
##' ## Use the 'autolayer' method for a quick comparison 
##' qM3 <-  quantMax(tv,
##'                  date = as.Date(sprintf("%4d-01-01", 2025:2035)),
##'                  level = c(0.95, 0.70))
##' gg + autolayer(qM3, colour = "SpringGreen3", linetype = "dashed") +
##'   ggtitle("Same as before. Green dashed line: restricted period 2025-2035")
##' ## Compare with a simulation. Note that 'sim' has class "bts"  and
##' ## is essentially a numeric matrix. Increase 'nsim' to get a more
##' ## precise estimate of the high quantiles
##' sim <- simulate(tv, newdate = date2, nsim = 10000)
##' M <- apply(sim, 2, max)
##' probs <- c(0.9, 0.95, 0.98, 0.99, 0.995, 0.998)
##' qSim <- quantile(M, prob = probs)
##' dfSim <- data.frame(Prob = probs, ProbExc = 1 - probs, Quant = qSim)
##' gg + geom_point(data = dfSim, mapping = aes(x = ProbExc, y = Quant))
##'\dontrun{
##' qM3 <- quantMax(tv, level = c(0.95, 0.70), confint = "proflik")
##' gg3 <- autoplot(qM3, fillConf = TRUE) +
##'                 ggtitle(paste("Quantile of the maximum over years 2025-2054 with",
##'                               " \"profile\" intervals"))
##' gg3
##' }
quantMax.TVGEV <- function(object,
                           prob,
                           date = NULL,
                           level = 0.95,
                           psi = NULL,
                           confintMethod = c("delta", "proflik"),
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
    } else if (any(prob <= 0) || any(prob >= 1.0)) {
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

        Quant <- array(NA_real_,
                       dim = c(Prob = nProb, Lim = 3L, Level = nLevel),
                       dimnames = list(Prob = fProb,
                                       Type = c("Quant", "L", "U"),
                                       Level = fLevel)) 
        
        ## =====================================================================
        ## Find a grid of quantiles which should cover the probability
        ## range.
        ## =====================================================================
   
        qMin <- .qMaxL.TVGEV(p = prob[1], theta)
        qMax <- .qMaxU.TVGEV(p = prob[length(prob)], theta)
        
        ## rq <- qMax - qMin
        ## qMin <- qMin - 0.1 * rq
        ## qMax <- qMax + 0.1 * rq
        mStarGrid <- seq(from = qMin, to = qMax, by = 0.001)

        ## BUG FIX: for a stationary model, 'qMin' and 'qMax' are
        ## equal and 'approx' would throw an error'.
        
        if (length(mStarGrid) > 2) {
            FMStarGrid <- rep(0.0, length(mStarGrid))
            
            for (i in seq_along(mStarGrid)) {
                FMStarGrid[i]<- prod(nieve::pGEV(mStarGrid[i],
                                                 loc = theta[ , 1],
                                                 scale = theta[ , 2],
                                                 shape = theta[ , 3]))
            }
            
            ## =================================================================
            ## An interpolation is required to get the values at
            ## the 'pretty' probability values.
            ## =================================================================
        
            qMStar <- approx(x = FMStarGrid, y = mStarGrid, xout = prob)$y
            
            if (DEBUG) {
                plot(prob, qMStar, pch = 16)
                lines(FMStarGrid, mStarGrid)
            }
            
        } else {
            qMStar <- qMin
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

    } else if (confintMethod == "proflik") {
        Quant <- array(NA_real_,
                       dim = c(Prob = nProb, Lim = 3L, Level = nLevel),
                       dimnames = list(Prob = fProb,
                                       Type = c("Quant", "L", "U"),
                                       Level = fLevel)) 
        diagno <- array(NA_real_,
                        dim = c(Prob = nProb, Type = 2L, Level = nLevel, Diag = 4L),
                        dimnames = list(Prob = fProb,
                                        Type = c("L", "U"),
                                        Level = fLevel,
                                        Diag = c("status", "objective", "constraint", "gradDist")))
        Psi <- array(NA_real_,
                     dim = c(Prob = nProb, Type = 2L, Level = nLevel, Param = p),
                     dimnames = list(Prob = fProb,
                                     Type = c("L", "U"),
                                     Level = fLevel,
                                     Param = parNames(object)))
        
        for (iProb in seq_along(prob)) {
            qFuni <- function(psi, object) {
                qMax.TVGEV(object = object, p = prob[iProb], date = date, psi = psi, deriv = TRUE,
                           trace = 1)
            }   
            resi <- profLik(object = object, fun = qFuni, level = level)
            
            Quant[iProb, , ] <- resi[c("est", "L", "U"), ]
            diagno[iProb, , , ] <- attr(resi, "diagno")
            Psi[iProb, , , ] <- attr(resi, "psi")
        }
        if (trace > 1) {
            cat("Diagnostics for the constrained optim\n")
            print(diagno)
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
        
    if (confintMethod == "proflik") attr(Quant, "psi") <- Psi
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
.qMaxL.TVGEV <- function(theta, p) {
    if (ncol(theta) != 3) {
        stop("'theta' must be a numeric matrix with three columns")
    }
    n <- nrow(theta)
    shapeMin <- min(theta[ , "shape"])
    scaleMin <- min(theta[ , "scale"])
    locMin <- min(theta[ , "loc"])
    nieve::qGEV(p = p^(1 / n),
                loc = locMin,
                scale = scaleMin,
                shape = shapeMin)
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
.qMaxU.TVGEV <- function(theta, p) {
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
##'     probability of non-exceedance for arbitrary given quantiles.
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
##'     parameters. The columns are in the order \emph{location},
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
##' date <- as.Date(sprintf("%4d-01-01", 2025:2054))
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

##' Compute the distribution function of the maximum over a period of
##' time.
##'
##' @title Value of the Distribution Function of the Maximum over a
##'     Period of Time
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param q A vector of "quantile" values at which the distribution
##'     function is to be evaluated.
##'
##' @param date A vector that can be coerced to the class
##'     \code{"Date"} which defines the period of time.
##'
##' @param psi Optional vector of parameters. By default, the
##'     estimated vector of parameters \code{coef(object)} is used.
##'
##' @param deriv Logical. If \code{TRUE} the gradient of the
##'     probability distribution will be computed and returned as the
##'     \code{"gradient"} attribute of the result.
##' 
##' @param hessian Logical. If \code{TRUE} the Hessian of the density
##'     will be computed and returned as the \code{"hessian"}
##'     attribute of the result.
##' 
##' @param ... Not used yet.
##'
##' @return A vector of values for the distribution function. Beside
##'     the gradient (w.r.t. the parameters) returned as the attribute
##'     named \code{"gradient"}, the Hessian w.r.t. \code{x} can
##'     also be returned as the attribute \code{"hessian"} of the result.
##' 
##' @export
##' 
pMax.TVGEV <- function(object,
                       q,
                       date = NULL,
                       psi = NULL,
                       deriv = FALSE,
                       hessian = FALSE,
                       ...) {

    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'deriv' is TRUE") 
    }
    
    if (is.null(psi)) psi <- object$estimate
    if (is.null(date)) date <- object$fDate
    date <- as.Date(date)
    theta <- psi2theta(model = object, psi = psi, date = date,
                       checkNames = FALSE)
    
    res <- rep(0.0, length(q))
    
    if (deriv) {
        mM <- modelMatrices(object, date = date)$X
        grad <- array(0.0, dim = c(length(q), object$p),
                      dimnames = list(paste0("q = ", format(q)), parNames(object)))
        if (hessian) {
            hess <- array(0.0, dim = c(length(q), object$p, object$p),
                          dimnames = list(paste0("q = ", format(q)),
                                          parNames(object),
                                          parNames(object)))
        }
    }
    
    for (i_q in seq_along(q)) {
        
        F_q <- nieve::pGEV(q[i_q],
                           loc = theta[ , 1],
                           scale = theta[ , 2],
                           shape = theta[ , 3],
                           deriv = deriv, hessian = hessian)
        
        res[i_q] <- prod(F_q)
        
        if (deriv) {
            
            ## 'grad_q'   is an array with dim c(length(date), 3)
            ## 'mM[[k]]'  is a matrix with dim c(length(date), p_k) where
            ##            'p_k' is the number of "psi" parameters corresponding to
            ##            'theta_k'
            grad_q <- attr(F_q, "gradient")
            grad_q <- sweep(grad_q, MARGIN = 1, STATS = F_q, FUN = "/")

            ## 'hess_q'   is an array with dimension (length(date), 3, 3)
            ## do the same for the Hessian H[b, k, ell] <- H[b, k, ell] / F[b]
            if (hessian) {
                hess_q <- attr(F_q, "hessian")
                hess_q <- sweep(hess_q, MARGIN = 1, STATS = F_q, FUN = "/")
            }

            ## Fill the gradient matrix
            for (i in 1:3) {
                ind_i <- object$ind[[i]]
                if (!object$isCst[i]) {
                    grad[i_q, ind_i] <- res[i_q] * (t(grad_q[ , i, drop = FALSE]) %*%
                                                    mM[[i]])
                } else {
                    grad[i_q, ind_i] <- res[i_q] * sum(grad_q[ , i])
                }
            }
            ## Fill the Hessian array
            if (hessian) {
                if (FALSE) {
                } else {
                    for (i in 1:3) {
                        ind_i <- object$ind[[i]]
                        for (j in 1:3) {
                            ind_j <- object$ind[[j]]
                            hm_ij <- 0.0
                            ## to be improved
                            for (b in 1:length(date)) {
                                hm_ij <- hm_ij +
                                    (hess_q[b, i, j] - grad_q[b, i] * grad_q[b, j]) *
                                    tcrossprod(mM[[i]][b, ], mM[[j]][b, ])
                            }
                            
                            A <- tcrossprod(grad[i_q, ind_i] /
                                            res[i_q], grad[i_q, ind_j] / res[i_q])
                            
                            A <- A + hm_ij
                            
                            hess[i_q, ind_i, ind_j] <- res[i_q] * A
                            
                            ## ## XXX
                            ## hess[i_q, ind_i, ind_j] <- hess[i_q, ind_i, ind_j] +
                            ##     tcrossprod(grad[i_q, ind_i] /
                            ##               res[i_q], grad[i_q, ind_j] / res[i_q])
                        }
                    }
                }
            }
        }
    }
    
    if (deriv) { 
        attr(res, "gradient") <- grad
        if (hessian)  attr(res, "hessian") <- hess
    }
    
    res
}

##' Compute the probability density of the maximum over a period of
##' time.
##'
##' @title Value of the Probability Density of the Maximum over a
##'     Period of Time
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param x A vector of values at which the density is to be
##'     evaluated.
##'
##' @param date A vector that can be coerced to the class
##'     \code{"Date"} which defines the period of time.
##'
##' @param psi Optional vector of parameters. By default, the
##'     estimated vector of parameters \code{coef(object)} is used.
##'
##' @param deriv Logical. If \code{TRUE} the gradient of the density
##'     will be computed and returned as the \code{"gradient"}
##'     attribute of the result.
##' 
##' @param derx Logical, If \code{TRUE} the derivative of the density
##'     w.r.t. \code{x} will be computed and returned as the
##'     \code{"derx"} attribute of the result.
##' 
##' @param ... Not used yet.
##'
##' @return A vector of values for the probability density. Beside the
##'     gradient (w.r.t. the parameters) returned as the attribute
##'     named \code{"gradient"}, the derivative w.r.t. \code{x} is
##'     also returned as the attribute \code{"derx"} of the result.
##' 
##' @export
##' 
dMax.TVGEV <- function(object,
                       x,
                       date = NULL,
                       psi = NULL,
                       deriv = FALSE,
                       derx = FALSE,
                       ...) {
    
    
    if (is.null(psi)) psi <- object$estimate
    if (is.null(date)) date <- object$fDate
    date <- as.Date(date)
    theta <- psi2theta(model = object, psi = psi, date = date,
                       checkNames = FALSE)
    
    if (deriv) {
        mM <- modelMatrices(object, date = date)$X
    }
    
    res <- rep(0.0, length(x))
    
    if (deriv) {
        grad <- array(0.0, dim = c(length(x), object$p),
                      dimnames = list(paste0("x = ", format(x)), parNames(object)))
        gradp <- array(0.0, dim = c(length(x), object$p),
                       dimnames = list(paste0("q = ", format(x)), parNames(object)))
    }
    
    if (derx) {
        gradx <- rep(0.0, length(x))
        names(gradx) <- paste0("x = ", format(x))
    }
    
    for (i_x in seq_along(x)) {
        
        FGEV_x <- nieve::pGEV(x[i_x],
                           loc = theta[ , 1],
                           scale = theta[ , 2],
                           shape = theta[ , 3],
                           deriv = deriv)
        
        fGEV_x <- nieve::dGEV(x[i_x],
                              loc = theta[ , 1],
                              scale = theta[ , 2],
                              shape = theta[ , 3],
                              deriv = deriv)
        
        FM  <- prod(FGEV_x)
        S1 <- sum(fGEV_x / FGEV_x)
        res[i_x] <- fM <- FM * S1
        
        if (deriv) {
            gradFGEV_x <- attr(FGEV_x, "gradient")
            gradFGEV_x <- sweep(gradFGEV_x, MARGIN = 1, STATS = FGEV_x, FUN = "/")
            
            mat <- sweep(attr(fGEV_x, "gradient"), MARGIN = 1,
                         STATS = FGEV_x, FUN = "/") -
                sweep(attr(FGEV_x, "gradient"), MARGIN = 1,
                      STATS = fGEV_x / FGEV_x / FGEV_x, FUN = "*")
            ## Then gradient w.r.t 'psi'
            grad0fM <- gradfM <- gradFM <- rep(NA_real_, object$p)
            
            for (i in 1:3) {
                ind_i <- object$ind[[i]]
                if (!object$isCst[i]) {
                    gradFM[ind_i] <- FM * (t(gradFGEV_x[ , i, drop = FALSE]) %*% mM[[i]])
                    grad0fM[ind_i] <- FM * t(mM[[i]]) %*% mat[ , i]
                } else {
                    gradFM[ind_i] <- FM * sum(gradFGEV_x[ , i])
                    grad0fM[ind_i] <- FM * sum(mat[ , i]) 
                }
            }
            
            ## prepare the second sum in dfM / dm
            gradfM <- gradFM * S1 + grad0fM
            grad[i_x, ] <- gradfM
            gradp[i_x, ] <- gradFM
        }

        if (derx) {
            ## prepare the second sum in dfM / dm
            S2  <- fGEV_x / FGEV_x            
            denom <-  theta[ , 2] + theta[ , 3] * (x[i_x] - theta[ , 1])
            ind <- denom > 0.0
            S2[ind] <- - S2[ind] * (theta[ind, 3] + 1.0) / denom[ind]
            S2 <- sum(S2[ind])
            gradx[i_x] <- fM * S1 + FM * S2 
        }
        
    }

    if (derx) attr(res, "derx") <- gradx

    if (deriv) {
        attr(res, "gradp") <- gradp
        attr(res, "gradient") <- grad
    }
    
    res
}

##' @title Quantile Function for a Random Maximum
##'
##' @param object An object that can be used to compute/predict the
##'     distribution of the random maximum.
##'
##' @param ... Further arguments for methods.
##' 
##' @return A function (or closure) that can be used to compute the
##'     quantiles for arbitrary given probabilities.
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
##' This function computes the quantile by using the distribution
##' function \eqn{F_{M^\star}}{F_M} returned by
##' \code{\link{cdfMaxFun.TVGEV}} and the \code{\link[stats]{uniroot}}
##' function to solve the equation \eqn{F_{M^\star}(q) = p}{F_M(q) =
##' p}. By contrast \code{quantMax.TVGEV} computes the values of the
##' distribution values \eqn{p_i} corresponding to a grid of quantiles
##' \eqn{q_i} and interpolates to find the quantiles corresponding to
##' the probability values provided by the user. So small differences
##' will exist in the results.
##' 
##' @title Quantile Function for the Random Maximum on a given
##'     Period
##' 
##' @param object An object with class \code{"TVGEV"}.
##' 
##' @param date An object that can be coerced to the class
##'     \code{"Date"}. If not provided, this will be taken as the date
##'     attached to \code{object}.
##'
##' @param psi An optional vector of coefficients for
##'     \code{object}. By default the ML estimate as returned by
##'     \code{coef(object)} is used.
##'
##' @param theta An optional matrix with three columns containing GEV
##'     parameters. The columns are in the order \emph{location},
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
##' @seealso \code{\link{cdfMaxFun.TVGEV}} for the corresponding
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
##' date <- as.Date(sprintf("%4d-01-01", 2025:2054))
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
            interv <- c(.qMaxL.TVGEV(theta = theta, p = probs[i]),
                        .qMaxU.TVGEV(theta = theta, p = probs[i]))
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

##' Compute the quantiles of the maximum over a period of time.
##'
##' @title Values of the Quantile Functon of the Maximum over a
##'     Period of Time
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param p A vector of probability values at which the quantile
##'     function is to be evaluated.
##'
##' @param date A vector that can be coerced to the class
##'     \code{"Date"} which defines the period of time.
##'
##' @param psi Optional vector of parameters. By default, the
##'     estimated vector of parameters \code{coef(object)} is used.
##'
##' @param deriv Logical. If \code{TRUE} the gradient of the quantile
##'     function will be computed and returned as the
##'     \code{"gradient"} attribute of the result.
##' 
##' @param hessian Logical, If \code{TRUE} the Hessian of the quantile
##'     function will be computed and returned as the \code{"hessian"}
##'     attribute of the result.
##'
##' @param trace Integer level of verbosity.
##' 
##' @param ... Not used yet.
##'
##' @return A vector of values for the quanitle function. Beside the
##'     gradient (w.r.t. the parameters) returned as the attribute
##'     named \code{"gradient"}, the Hessian can also be returned as
##'     the attribute \code{"Hessian"} of the result.
##' 
##' @export
##'
##' @note The \code{\link{quantMax}} method can be used to compute the
##'     quantiles and confidence intervals: delta or profile
##'     likelihood.
##'
##' @seealso The \code{\link{quantMax.TVGEV}} method.
##' 
qMax.TVGEV <- function(object,
                       p,
                       date = NULL,
                       psi = NULL,
                       deriv = FALSE,
                       hessian = FALSE,
                       trace = 0,
                       ...) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'deriv' is TRUE") 
    }

    if (is.null(psi)) psi <- object$estimate
    if (is.null(date)) date <- object$fDate
    date <- as.Date(date)
    theta <- psi2theta(model = object, psi = psi, date = date,
                       checkNames = FALSE)
    sigmaMin <- min(theta[ , 2])

    if (sigmaMin <= 0.0) {
        return(NA)
    }
        
    FMaxZero <- function(q, prob) {
        pMax.TVGEV(q, date = date, psi = psi, object = object) - p
    }
        
    res <- rep(0.0, length(q))
    
    if (deriv) {
        grad <- array(0.0, dim = c(length(p), object$p),
                      dimnames = list(paste0("p = ", format(p)),
                                      parNames(object)))
        if (hessian) {
            warning("The computation of the Hessian is experimental. ",
                    "Though the formulas used seem OK, some instability ",
                    "in the computation may arise when p >= 0.9")
            hess <- array(0.0, dim = c(length(p), object$p, object$p),
                          dimnames = list(paste0("p = ", format(p)),
                                          parNames(object),
                                          parNames(object)))
            mM <- modelMatrices(object, date = date)$X
        } else {
            hess <- NULL
        }
    } else {
        grad <- NULL
    }
    
    for (i_p in seq_along(p)) {
        
        interv <- c(.qMaxL.TVGEV(theta = theta, p = p[i_p]),
                    .qMaxU.TVGEV(theta = theta, p = p[i_p]))

        ## An error would result when the interval is one point.  This
        ## will happen when 'date' has length one or when the model is
        ## stationary.
        ##
        ## XXX Add such models to the tests.
        
        cond <- (interv[1] == interv[2])
        
        if (interv[1] == interv[2]) {
            res[i_p] <- interv[1]
        } else {
            resTry <- try(uniroot(f = FMaxZero,  interval = interv,
                                  prob = p[i_p]))
            ## res[i] is the quantile
            if (!inherits(resTry, "try-error")) {
            res[i_p] <- resTry$root
            }
        }
        
        if (deriv) {

            ## ===============================================================
            ## we need to re-compute the values of the GEV density and
            ## of the GEV distribution function since these are not
            ## stored. Maybe some savings could later be obtained by
            ## optionnaly storing these?
            ## ===============================================================
            
            fGEV_p <- nieve::dGEV(res[i_p],
                               loc = theta[ , 1],
                               scale = theta[ , 2],
                               shape = theta[ , 3], deriv = hessian)
            FGEV_p <- nieve::pGEV(res[i_p],
                               loc = theta[ , 1],
                               scale = theta[ , 2],
                               shape = theta[ , 3], deriv = hessian)
            
            FM <- pMax.TVGEV(res[i_p], date = date, psi = psi, object = object,
                             deriv = TRUE, hessian = hessian)

            sumLog <- sum(fGEV_p / FGEV_p)
            fM <- FM * sumLog
            
            grad[i_p, ] <- -attr(FM, "gradient") / fM

            if (hessian) {

                ## =============================================================
                ## 'grad0fM' is the gradient of the density computed
                ## by ignoring the dependence of 'm' on 'psi' in the
                ## differentiation. The computation here is similar to
                ## that done in the 'dMax.TVGEV' function. The same is
                ## true for 'derfM' which is coputed as in
                ## 'dMax.TVGEV' when 'derx' is TRUE.
                ## =============================================================
                
                ## we need to compute the gradient of  
                ## First gradient w.r.t. the GEV parameter 'theta'.
                ## 'mat' is a Bstar * 3 matrix
                mat <- sweep(attr(fGEV_p, "gradient"), MARGIN = 1,
                             STATS = FGEV_p, FUN = "/") -
                    sweep(attr(FGEV_p, "gradient"), MARGIN = 1,
                          STATS = fGEV_p / FGEV_p / FGEV_p, FUN = "*")
                
                ## Then gradient w.r.t 'psi'
                grad0fM <- rep(NA_real_, object$p)
                
                for (i in 1:3) {
                    ind_i <- object$ind[[i]]
                    if (!object$isCst[i]) {
                        grad0fM[ind_i]  <- FM * t(mM[[i]]) %*% mat[ , i]
                    } else {
                        grad0fM[ind_i] <- FM * sum(mat[ , i]) 
                    }
                }

                gradFM <- attr(FM, "gradient") ## to get an easier-to-read code
                grad0fM <- gradFM * sumLog + grad0fM 

                ## 'derLogfGEV_p' is a vector with length Bstar := nrow(theta)
                derLogfGEV_p <- fGEV_p / FGEV_p
                denom <-  theta[ , 2] + theta[ , 3] * (res[i_p] - theta[ , 1])
                ind <- denom > 0.0
                derLogfGEV_p[ind] <- - derLogfGEV_p[ind] *
                    (theta[ , 3] + 1.0) / denom
                derfM <- fM * sumLog + FM * sum(derLogfGEV_p)

                ## This was wrong
                ## ## Correction of the gradient because 'm' depends on 'psi'
                ## if (trace) {
                ##     cat("correction of grad 'f_M' due to the dependence of",
                ##         " 'm' on 'psi'\n")
                ##     print(- derfM * gradFM / fM)
                ## }
                ## gradfM <- - derfM * gradFM / fM + grad0fM
                
                gradfM <- grad0fM
                
                if (trace > 1) {
                    cat("\nValue of F_M\n")
                    FMP <- FM
                    attributes(FMP) <- NULL
                    print(FMP)
                    
                    cat("\nGradient of f_M\n")
                    print(gradfM)
                    
                    cat("\nGradient of F_M\n")
                    print(gradFM)
                   
                    fMCheck <- dMax.TVGEV(object, date = date,
                                          x = res[i_p], deriv = TRUE, derx = TRUE)
                    cat("\nCheck the density 'f_M'\n")
                    print(c("here" = as.numeric(fM),
                            "check" = as.numeric(fMCheck)))
                    cat("\nCheck the derivative 'derf_M'\n")
                    print(c("here" = as.numeric(derfM),
                            "check" = as.numeric(attr(fMCheck, "derx"))))
                }

                ## =============================================================
                ## compute the derivatives of fM w.r.t the 1-st arg and
                ## w.r.t. the parameters
                ## Mind that both 'gradfM' and 'gradFM' are row matrices
                ## so 'crossprod' acts as a 'tcrossprod' on vectors.
                ## =============================================================
                
                hess[i_p, , ] <- (-attr(FM, "hessian")[1, , ] +
                                  (crossprod(gradfM, gradFM) +
                                   crossprod(gradFM, gradfM) -
                                   derfM * crossprod(gradFM) / fM) / fM) / fM

                if (trace > 1) {
                    cat("\nComponents in Hessian\n")
                    print(-attr(FM, "hessian")[1, , ] / fM)
                    print((crossprod(gradfM, gradFM) +
                           crossprod(gradFM, gradfM)) / fM / fM)
                    print(-derfM * crossprod(gradFM) / fM / fM / fM)       
                }
                
            }   
        }
    }
    
    if (deriv) {
        attr(res, "gradient") <- grad
        if (hessian) {
            attr(res, "hessian") <- hess
        }
    }
    
    res
    
}

    
