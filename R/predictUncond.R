## ************************************************************************
##' Prediction of unconditional or Non-Stationary Return Levels for a
##' \code{NSGEV} or a \code{TVGEV} object.
##'
##' For a Non-Stationary GEV model, the Return Level corresponding to
##' a given period (e.g. 100 years) depends on the value of the
##' covariates. However we can make an assumption on the the
##' distribution of the covariates to compute an uncondtitional Return
##' Level as the expectation of the RL. Similarily, for Time-Varying
##' \code{TVGEV} models we can define the Return Level as the level
##' for which the expected number of exceedances on a reference period
##' of time is equal to one, see \code{\link{RL}}. In both cases, the RL
##' no longer depends of the covariates or on a specific block.
##' 
##' @title Prediction of Unconditional or Non-Stationary Return Levels
##' for a \code{NSGEV} object or \code{TVGEV} object
##'
##' @param object A \code{NSGEV} or \code{TVGEV} model.
##'
##' @param period A numeric vector giving the periods at which the
##' Return Levels will be computed.
##'
##' @param level The Confidence Levels at which the CI will be
##' computed.  For \code{"proflik"} CI, only one level can be given
##' for now.
##'
##' @param newdata Only for a \code{NSGEV} model. A data frame
##' containing the covariates needed. Each row represents a block with
##' unit duration.
##' 
##' @param newdate Only for a \code{TVGEV} model. A vector that can
##' be coerced to the \code{"Date"} class.
##'
##' @param newdateFrom Only for a \code{TVGEV} model when
##' \code{newdate} is not given. A vector with length one that can be
##' coerced to the \code{"Date"} class. This gives the first block of a period
##' of \eqn{B} successive blocks used to compute the Return Levels,
##' where \eqn{B} is the maximal period defined in \code{period}.
##'
##' @param resample Not used yet. See the \code{\link{RL}}
##' function. For now, the RL period that are greater than
##' \code{nrow(newdata)} are discarded in the computation.
##'
##' @param RLType The type of Return Level as in \code{\link{RL}}.
##'
##' @param confintMethod The type of Confidence Interval (CI) to
##' compute.  The value \code{"delta"} leads to using the \emph{delta
##' method}, and \code{"proflik"} leads to using the
##' profile-likelihood method.
##'
##' @param out Type of output.
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
##' @section Caution: When \code{period} contains some large values,
##' say 100 or more, a massive extrapolation from the model is
##' required to compute the correspondint RLs. This is likely to be
##' meaningless from a statistical point of view, because a massive
##' extrapolation of regression-like model does not make sense. The
##' interest should rather focus on a moderately large period.
##' 
predictUncond <- function(object,
                          period = NULL,
                          level = 0.95,
                          newdata = NULL,
                          newdate = NULL,
                          newdateFrom = NULL,
                          resample = FALSE,
                          RLType = c("exceed", "average"), 
                          confintMethod = c("delta", "none", "proflik"),
                          out = c("data.frame", "array"),
                          trace = 0,
                          ...) {
    
    RLType <- match.arg(RLType)
    confintMethod <- match.arg(confintMethod)
    out <- match.arg(out)
    
    if (is.null(period)) period <- c(5, 10, 20, 30, 40, 50, 70, 100)
    else period <- sort(period)
    np <- length(period)
    
    if (is(object, "NSGEV")) {
        
        if (!is.null(newdate) || !is.null(newdateFrom)) {
            stop("'newdate' and 'newdateFrom' can only be used\n",
                 "with a TVGEV object, not with a NSGEV.")
        }

        if (is.null(newdata)) newdata <- object$data
        nd <- nrow(newdata)
        if (!resample && period[np] > nd) {
            warning("When 'resample' is FALSE, no prediction can be given ",
                    "for return periods greater than nrow(newdata)")
        }
    } else if (is(object, "TVGEV")) {
        
        if (!is.null(newdata)) {
            stop("'newdata' can only be used\n",
                 "with a NSGEV object, not with a TVGEV.")
        }
        
        if (is.null(newdate)) {
            if (is.null(newdateFrom )) {
                newdate <- as.Date(object$fDate)
            } else {
                newdate <- seq(from = as.Date(newdateFrom ),
                               by = "year", length.out = period[np])
            }
            nd <- length(newdate)
        }
    }

    if (confintMethod == "none") {

        RL <- array(NA, dim = c(np, 1),
                    dimnames = list("period" = period,
                        "Quant"))
                       
        for (i in seq_along(period)) {
            if (period[i] <= nd) {
                res <- RL(model = object,
                          period = period[i],
                          data = newdata,
                          date = newdate,
                          type = RLType, deriv = FALSE)
                RL[i, "Quant"] <- res
            }
        }
        
    } else {
        
        if (any(level >= 1.0)) stop("values in 'level' must be < 1.0")
        level <- sort(level)
   
        RL <- array(NA, dim = c(np, 3, length(level)),
                    dimnames = list("period" = period,
                        c("Quant", "L", "U"),
                        level))
    }
    
    if (confintMethod == "delta") {

        if (is.null(object$vcov)) {
            stop("'object' must have a 'vcov' element")
        }

        alpha <- 1 - level
        probsCI <- rep(NA, 2 * length(level))
        probsCI[seq(from = 1L, to = 2 * length(level), by = 2)] <- alpha / 2  
        probsCI[seq(from = 2L, to = 2 * length(level), by = 2)] <-  1 - alpha / 2
        qCI <- qnorm(probsCI, mean = 0, sd = 1)
        
        for (i in seq_along(period)) {
            if (period[i] <= nd) {
                res <- RL(model = object,
                          period = period[i],
                          data = newdata,
                          date = newdate,
                          type = RLType,
                          deriv = TRUE)
                grad <- attr(res, "gradient")
                RL[i, "Quant", 1L] <- res
                sd <- as.vector(grad %*% object$vcov %*% t(grad))
                RL[i , c("L", "U"), 1L] <- res + qCI * sd
            }
        }
        
    } else if (confintMethod == "proflik") {
        
        if (length(level) > 1L) {
            stop("when 'confInt' is \"proflik\", 'confLev' must be ",
                 " of length one (for now)")
        }   
        
        for (iPer in seq_along(period)) {
            if (period[iPer] <= nd) {

                if (trace) {
                    cat(sprintf("\n   o period %d\n   ============\n",
                                period[iPer]))
                }
                
                RLfun <- function(psi, object) {
                    rl <- RL(period = period[iPer], model = object,
                             data = newdata,
                             date = newdate,
                             psi = psi,
                             type = RLType,
                             deriv = TRUE)
                    if (is.na(rl)) {
                        rl <- NaN
                        attr(rl, "gradient") <- rep(NaN, object$p)
                    }
                    rl 
                }
                pl <- try(profLik(object = object, fun = RLfun, level = level, deriv = TRUE,
                                  trace = trace))
                if (!inherits(pl, "try-error")) {
                    RL[iPer, , 1L] <- pl[ , 1L]
                }
             }   
        }

    }
    
    if (out == "data.frame") {
        
        if (requireNamespace("reshape2", quietly = TRUE)) {
            
            if (confintMethod == "none") {
                RL <- data.frame(Period = period, RL)
            } else {
                ## ====================================================================
                ## UGGLY CODE: there must be a simpler and more
                ## efficient way of doing that. The problem is that we
                ## want to drop the "Type" dimension but not the
                ## "Level" dimension even when it has no extension
                ## ====================================================================
                
                df <- list()
                for (nm in c("Quant", "L", "U")) {
                    RL1 <- RL[ , nm, , drop = FALSE]
                    df[[nm]]  <- melt(RL1,
                                      value.name = nm,
                                      varnames = c("Period", "Type", "Level"))
                }
                RL <- data.frame(df[["Quant"]][, c("Period", "Level", "Quant")],
                                 L = df[["L"]]$L, U = df[["U"]]$U)
                RL$Level <- as.factor(RL$Level)
                levels(RL$Level) <- formatPerc(level)
                ## RL$Date <- as.Date(RL$Date)
                ind <- with(RL, order(Level, Period))
                RL <- RL[ind, , drop = FALSE]
            }
        } else {
            stop("the package 'reshape2' could not be used")
        }
        
        if (is(object, "TVGEV")) {
            RL <- data.frame(Date = newdate[1L], RL)
            class(RL) <- c("predict.TVGEV", "data.frame")
        }
    }
    
    return(RL)
    
}
