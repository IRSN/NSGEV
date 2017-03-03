##' Prediction of Return Levels for a \code{NSGEV} object
##'
##' @title Prediction of Return Levels for a \code{NSGEV} object
##'
##' @param object A \code{NSGEV} model.
##'
##' @param period A numeric vector giving the periods at which the
##' Return Levels will be computed.
##'
##' @param newdata A data frame containing the covariates needed. Each
##' row represents a block with unit duration.
##'
##' @param resample Not used yet. See the \code{\link{RL}}
##' function. For now, the RL period that are greater than
##' \code{nrow(newdata)} are discarded in the computation.
##'
##' @param RLType The type of Return Level as in \code{\link{RL}}.
##'
##' @param confInt The type of Confidence Interval (CI) to compute.
##' The value \code{"delta"} leads to using the \emph{delta method},
##' and \code{"proflik"} leads to using the profile-likelihood method.
##' 
##' @param confLevel The Confidence Levels at which the CI will be
##' computed.  For \code{"proflik"} CI, only one level can be given
##' for now.
##'
##' @param trace Level of verbosity.
##'
##' @param ... Not used yet.
##'
##' @return A list with several elements.
##' 
##' \item{RL}{A data frame with the predicted RLs and the related CIs,
##' with one row by block.}
##' 
##' \item{PsiStar}{When \code{confInt} is equal to \code{"proflik"},
##' this element is a matrix with its row \eqn{i} giving the value
##' \eqn{\boldsymbol{psi}^\star}{\psi*} of the vector of parameters
##' that maximizes the Return Level \eqn{\rho(T)} with period \eqn{T =
##' T_i} under the constraint on the log-likelihood. Since most often
##' the rows are close enough, a significant reduction of the
##' computing time could be achieved in the near future by using the
##' same value of \eqn{\boldsymbol{\psi}}{\psi} for all the Return
##' Periods.  }
##'  
##' @author Yves Deville
##'
##' @note For the profile-likelihood method, the determination of the
##' confidence intervals is quite slow because a constrained
##' optimization problem is solved for each period.
##'
##' Future versions might allow different durations across blocks by
##' using dedicated arguments. For now it must be kept in mind that
##' the periods are understood as multiple of a constant block
##' duration. So if \code{newdata} has \code{100} rows the maximal
##' Return Period that can be used without resampling is \code{100}.
##'
##' @examples
##' example(as.NSGEV.fevd)
##' L.delta <- predict(ns, period = seq(from = 10,  to = 68, by = 10),
##'              RLType = "exceed")
##' plot(U_95 ~ period, data = L.delta$RL, ylim = c(20, 27), type = "l")
##' lines(quant ~ period, data = L.delta$RL)
##' 
##' ## much slower, but much trustier!
##' \dontrun{
##' L.proflik <- predict(ns, period = seq(from = 10,  to = 68, by = 10),
##'              RLType = "exceed", confInt = "proflik") 
##' plot(U_95 ~ period, data = L.proflik$RL, ylim = c(20, 30), type = "l")
##' lines(quant ~ period, data = L.proflik$RL)
##' }
predict.NSGEV <- function(object, period = NULL,
                          newdata = NULL, 
                          resample = FALSE,
                          RLType = c("exceed", "average"),
                          confInt = c("delta", "none", "proflik"),
                          confLevel = 0.95,
                          trace = 0,
                          ...) {

    RLType <- match.arg(RLType)
    confInt <- match.arg(confInt)
    
    if (is.null(newdata)) newdata <- object$data
    if (is.null(object$response)) {
        stop("the NSGEV model given in 'object' must have a reponse")
    } else {
        y <- object$response
    }

    if (is.null(period)) {
        period <- c(5, 10, 20, 50, 100, 150, 200, 500, 1000)
    } else {
        period <- sort(period)
        if (!resample && any(period > nrow(newdata))) {
            warning("When 'resample' is FALSE, no prediction can be given ",
                    "for return periods greater than nrow(newdata)")
        }
    }

    np <- length(period)
    pred <- data.frame(period = period, quant = rep(NA, np))

    if (confInt == "none") {
        for (i in seq_along(period)) {
            if (period[i] <= nrow(newdata)) {
                res <- RL(model = object, period = period[i], data = newdata,
                          type = RLType, deriv = FALSE)
                pred[i, "quant"] <- res
            }
        }

        return(list(RL = pred))
        
    } else {

        ##======================================================================
        ## Prepare a data frame with suitable names
        ##======================================================================
        
        if (any(confLevel >= 1.0)) stop("values in 'confLevel' must be < 1.0")
        confLevel <- sort(confLevel)
        confL <- 100 * confLevel
        dL <- confL - floor(confL)
        confNames <- sprintf("_%2.0f", confL)
        ind <- (dL > 1e-3)
        if (any(ind))  confNames[ind] <- sprintf("_%3.1f", confL[ind])
        confNames <- paste(rep(c("L", "U"), times = length(confNames)),
                           rep(confNames, each = 2L), sep = "")
        alpha <- 1 - confLevel
        probsCI <- rep(NA, 2 * length(confLevel))
        probsCI[seq(from = 1, to = 2 * length(confLevel), by = 2)] <- alpha / 2  
        probsCI[seq(from = 2, to = 2 * length(confLevel), by = 2)] <-  1 - alpha / 2
        qCI <- qnorm(probsCI, mean = 0, sd = 1)
        CL <- array(NA, dim = c(np, 2 * length(confLevel)),
                    dimnames = list(period, confNames)) 
    }
    
    if (confInt == "delta") {

        if (is.null(object$vcov)) {
            stop("'object' must have a 'vcov' element")
        }
        
        for (i in seq_along(period)) {
            if (period[i] <= nrow(newdata)) {
                res <- RL(model = object, period = period[i], data = newdata,
                          type = RLType, deriv = TRUE)
                grad <- attr(res, "gradient")
                pred[i, "quant"] <- res
                sd <- grad %*% object$vcov %*% t(grad)
                CL[i , ] <- res + qCI * sd
            }
        }

        pred <- cbind(pred, CL)
        return(list(RL = pred))
        
    } else if (confInt == "proflik") {

        
        if (length(confLevel) > 1L) {
            stop("when 'confInt' is \"proflik\", 'confLev' must be ",
                 " of length one (for now)")
        }   
        
        ## this is the NEGATIVE LogLik
        ellL  <-  object$negLogLik + qchisq(confLevel, df = 1) / 2.0
        
        ##=====================================================================
        ## This is the objective function, shipping the gradient with
        ## the result.  The opposite of the RL since we are
        ## minimizing!
        ## ====================================================================
        f <- function(psi, period) {
            res <- RL(period, model = object, data = newdata, psi = psi, deriv = TRUE)
            list("objective" = -res, "gradient" = -attr(res, "gradient"))
        }
        
        ##=====================================================================
        ## This is the constraint function, shipping the gradient with
        ## the result.  -ell(psi) + ell_max - delta
        ## ====================================================================
        g <- function(psi, period) {
            res <- negLogLik(psi = psi, model = object, data = newdata,
                             deriv = TRUE, checkNames = FALSE)
            theta <- psi2theta(model = object, psi = psi, data = newdata,
                               deriv = TRUE, checkNames = FALSE)
            list("constraints" = res$objective - ellL, "jacobian" = res$gradient)
        }
        
        opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel" = 1.0e-10,
                      "maxeval" = 1000,
                      ## "check_derivatives" = TRUE, "check_derivatives_print" = "all",
                      "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7),
                      "print_level" = 0)
        
        ## structures to save the results
        res <- list()
        PsiStar <- array(NA, dim = c(length(period), object$p),
                         dimnames = list(period, parNames(object)))
     
        Rho <- rep(NA, length(period))
        names(Rho) <- period
        
        ## Initial 
        psi <- object$estimates
        
        for (i in seq_along(period)) {
            if (period[i] <= nrow(newdata)) {
                rl <- RL(model = object, period = period[i], data = newdata,
                          type = RLType, deriv = FALSE)
                
                res[[i]] <- nloptr(x0 = psi,
                                   eval_f = f,
                                   eval_g_ineq = g,
                                   period = period[i],
                                   opts = opts1)
                
                psi <- res[[i]][["solution"]]
                Rho[i] <- -res[[i]][["objective"]]
                PsiStar[i, ] <- psi
                pred[i, "quant"] <- rl
                CL[i, 2L] <- Rho[i]
             }   
        }
        
        pred <- cbind(pred, CL)
        print(PsiStar)
        return(list(RL = pred, PsiStar = PsiStar))

    }

}
