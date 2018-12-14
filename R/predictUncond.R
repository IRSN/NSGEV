## ************************************************************************
##' Prediction of Unconditional or Non-Stationary Return Levels for a
##' \code{TVGEV} or a \code{NSGEV} object.
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
##' for a \code{TVGEV} or a \code{NSGEV} Object
##'
##' @param object A \code{TVGEV} or a \code{NSGEV} model object.
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
##' @param newdate Only for a \code{TVGEV} model. This argument is
##' essentially for compatibility with the \code{predict} method for
##' unconditional prediction. A vector that can be coerced to the
##' \code{"Date"} class. This vector define \eqn{B} blocks used to
##' compute the return levels. The standard use is is with consecutive
##' blocks, so the specification should \emph{use the}
##' \code{newdateFrom} \emph{formal argument rather than} \code{newdate}.
##'
##' @param newdateFrom Only for a \code{TVGEV} model. A vector with
##' length one that can be coerced to the \code{"Date"} class. This
##' gives the first block of a period of \eqn{B} successive blocks
##' used to compute the Return Levels, where \eqn{B} is the maximal
##' period defined in \code{period}.
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
##' @return A data frame with the predicted RLs and the related CIs,
##' with one row by block.
##'
##' This data frame can have several attributes
##'
##' \itemize{
##'
##' \item \code{title} is a string that can be used as title in plots.
##'
##' \item \code{diagno}. When \code{confInt} is equal to
##' \code{"proflik"}, this element is an array providing information
##' about the optimisation.
##' 
##' \item \code{PsiStar}. When \code{confInt} is equal to
##' \code{"proflik"}, this element is a matrix with its row \eqn{i}
##' giving the value \eqn{\boldsymbol{psi}^\star}{\psi*} of the vector
##' of parameters that maximizes the Return Level \eqn{\rho(T)} with
##' period \eqn{T = T_i} under the constraint on the log-likelihood.
##'
##' }
##'  
##' @author Yves Deville
##'
##' @note For the profile-likelihood method, the determination of the
##' confidence intervals is quite slow because a constrained
##' optimization problem is solved for each period.
##'
##' Since most often the values of \eqn{\boldsymbol{\psi}^\star}{\psi*}
##' corresponding to different periods are close enough, a significant
##' reduction of the computing time could be achieved in the near
##' future by using the same value of \eqn{\boldsymbol{\psi}}{\psi}
##' for all the Return Periods.
##'
##' Future versions might allow different durations across blocks by
##' using dedicated arguments. For now, it must be kept in mind that
##' the periods are understood as multiple of a constant block
##' duration. So if \code{newdata} has \code{100} rows the maximal
##' Return Period that can be used without resampling is \code{100}.
##' 
##' @section Caution:
##' \itemize{
##'
##' \item When \code{period} contains some large values, say 100 or
##' more, a massive extrapolation from the model is required to
##' compute the corresponding RLs. This is likely to be meaningless
##' from a statistical point of view, because a massive extrapolation
##' of regression-like model does not make sense. The interest should
##' rather focus on a moderately large period.
##'
##' \item For each period \eqn{T}, the Return Level is computed using
##' the \eqn{T} blocks \eqn{t_0}, \eqn{t_0 + 1}, \dots, \eqn{t_0 + T -1}
##' so the different periods \eqn{T} appearing on a Return Level plot
##' correspond to different prediction horizons. As a result, the
##' quantile or confidence bounds do not depend as smoothly of \eqn{T}
##' as they do for conditional predictions.
##'
##' }
##'
##' @examples
##' example(TVGEV)
##' pu <- predictUncond(res2, newdateFrom = "2020-01-01")
##' pu
predictUncond <- function(object,
                          period = NULL,
                          level = 0.95,
                          newdata = NULL,
                          newdate = NULL,
                          newdateFrom = NULL,
                          resample = FALSE,
                          RLType = c("exceed", "average"), 
                          confintMethod = c("delta", "none", "boot", "proflik"),
                          out = c("data.frame", "array"),
                          trace = 0,
                          ...) {

    
    probL <- (1 - level) / 2
    probU <- 1 - probL

    dots <- match.call(expand.dots = FALSE)$...
    mc <- match.call()
    ## cat("dots \n")
    ## print(dots)
    
    RLType <- match.arg(RLType)
    confintMethod <- match.arg(confintMethod)
    out <- match.arg(out)

    if (length(dots)) {
        dotsText <- paste(sprintf("'%s'", names(dots)), collapse = ", ")
        if (confintMethod %in% c("none", "delta", "proflik")) {
            warning("the formals ", dotsText, " are ignored")
        } else if (confintMethod == "boot") {
            fm <- formals(getS3method("bs", class(object)))
            ind <- !(names(dots) %in% names(fm))
            if (any(ind)) {
                dotsText2 <- paste(sprintf("'%s'", names(dots)[ind]),
                                   collapse = ", ")
                stop("the formals ", dotsText2, " are passed to 'bs'\n",
                     "method but are nor formals of it.")
            }
        }
    }
    
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
            
            if (is.null(newdateFrom)) {
                newdate <- as.Date(object$fDate)
                newdateFrom <- newdate[1] 
            } else {
                newdate <- seq(from = as.Date(newdateFrom),
                               by = "year", length.out = period[np])
            }
            
        } else {
            warning("For unconditional predictions, the normal use is via ",
                    "the formal argument 'newdateFrom' to specify consecutive ",
                    "blocks. Using the 'newdate' argument can be misleading.")
                    
        }
        nd <- length(newdate)
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

        diagno <- Psi <- NULL
        
    } else {

        overPeriod <- FALSE
        
        if (trace) {

            cat(sprintf(paste("The confidence bound are obtained by the '%s'",
                              " method.\n"), confintMethod))
        }
        
        if (any(level >= 1.0)) stop("values in 'level' must be < 1.0")
        level <- sort(level)
   
        RL <- array(NA, dim = c(np, 3, length(level)),
                    dimnames = list("period" = period,
                        c("Quant", "L", "U"),
                        level))

        for (i in seq_along(period)) {
            if (period[i] <= nd) {
                res <- RL(model = object,
                          period = period[i],
                          data = newdata,
                          date = newdate,
                          type = RLType, deriv = FALSE)
                RL[i, "Quant", ] <- res
            }
        }
        
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
                sd <- sqrt(as.vector(grad %*% object$vcov %*% t(grad)))
                RL[i , c("L", "U"),  ] <- res + qCI * sd
            }
        }

        diagno <- Psi <- NULL
        
    } else if (confintMethod == "boot") {
        
        if (is.null(object$boot)) {
            if (trace) {
                cat("Since 'object' does not embed a bootstrap distribution, \n",
                    "we use `bs` to compute this by bootstraping model\n",
                    "parameters")
                
            }
            object$boot <- bs(object, ...)
            
        } else {
            if (length(dots)) {
                message("'object' embeds a bootstrap distribution. Formal arguments\n",
                        " intended to be passed to `bs` will be ignored")
            }
        }
       
        RLloc <- function(psi, period) {
            RL(period = period, model = object,
               data = newdata, date = newdate, psi = psi,
               type = RLType, deriv = FALSE)
        }

        for (iPer in seq_along(period)) {
            if (period[iPer] <= nd) {
                RLs <- apply(X = object$boot$estimate, MARGIN = 1L,
                             FUN = RLloc, period = period[iPer])
                RL[iPer, "L",  ] <- quantile(RLs, prob = probL)
                RL[iPer, "U",  ] <- quantile(RLs, prob = probU)
            } else if (!overPeriod) {
                overPeriod <- TRUE
            }
                
        }
        
        diagno <- Psi <- NULL
        
    } else if (confintMethod == "proflik") {
        
        ## if (length(level) > 1L) {
        ##     stop("when 'confInt' is \"proflik\", 'confLev' must be ",
        ##          " of length one (for now)")
        ## }
        
        diagno <-
            array(NA,
                  dim = c(length(period), 2L, length(level), 4L),
                  dimnames = list(period, c("L", "U"), level,
                      c("status", "objective", "constraint", "gradDist")))
        Psi <-
            array(NA,
                  dim = c(length(period), 2L, length(level),  object$p),
                  dimnames = list(period, c("L", "U"), level, object$parNames))
        
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
                
                pl <- try(profLik(object = object, fun = RLfun,
                                  level = level, deriv = TRUE,
                                  trace = trace))
                
                if (!inherits(pl, "try-error")) {
                    RL[iPer, , ] <- pl[ , ]
                    diagno[iPer, , , ] <- attr(pl, "diagno")
                    Psi[iPer, , , ] <- attr(pl, "psi")
                }
            }   
        }
       
    }
    
    if (out == "data.frame") {
        
        if (requireNamespace("reshape2", quietly = TRUE)) {
            
            if (confintMethod == "none") {
                RL <- data.frame(Period = period, RL)
            } else {
                ## ============================================================
                ## UGGLY CODE: there must be a simpler and more
                ## efficient way of doing that. The problem is that we
                ## want to drop the "Type" dimension but not the
                ## "Level" dimension even when it has no extension
                ## ============================================================
                
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
    attr(RL, "diagno") <- diagno
    attr(RL, "psi") <- Psi
    attr(RL, "title") <-
        sprintf("Integrated Return Levels for periods starting at  %s",
                format(newdateFrom))
    attr(RL, "type") <- "unconditional"
    invisible(RL)
    
}
