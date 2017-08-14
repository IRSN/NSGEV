## ************************************************************************
##' Predict method for \code{TVGEV} objects.
##'
##' Compute the Return Levels (RLs) for different periods and
##' different values of the time block, as well as Confidence
##' Intervals (CIs) for these. The results are not predictions in the
##' usual acceptance.  The CIs can be obtained by the usual
##' \emph{delta method}, provided that the approximate covariance for
##' the estimated parameters is found in \code{object}. They also can
##' be obtained by \emph{bootstrap}, using by default the bootstrap
##' distribution embeded in \code{object} if any or by computing it
##' else. The bootstrap can be parametric or non-parametric, see
##' \code{\link{bs.TVGEV}}. Finally, the \emph{profile-likelihood}
##' method can be used: the confidence limits for a RL are obtained by
##' maximising and minimising this RL under the constraint that the
##' log-likelihood is greater than a suitable value. This method
##' avoids the usual re-parameterisation of the model which can be
##' tedious for the general form of model allowed in \code{TVGEV}.
##' 
##' @title  Predict Method for \code{TVGEV} Objects
##'
##' @param object An object with S3 class \code{"TVGEV"} 
##'
##' @param newdate A vector with class \code{"Date"} or that can be
##' coerced to this class.
##'
##' @param period Numeric vector of periods expressed as multiple of
##' the block duration. Usually the block duration is one year and a
##' value \code{100} will correspond to a 100-year return level, i.e.
##' to the GEV quantile with probability 0.99. The default value
##' allows the construction of an acceptable RL plot, but if a
##' smoother RL curve and/or Confidence Band is needed, more values
##' would be used.
##'
##' @param level Numeric vector of confidence level(s).
##'
##' @param confintMethod The method used to compute the confidence
##' intervals. See \bold{Details}.
##'
##' @param out Type of output.
##' 
##' @param biasCorrect \code{Logical} used only when
##' \code{confintMethod} is \code{"boot"}. When \code{TRUE}, the RL
##' named \code{"Quantile"} is computed as the average of the RLs
##' obtained with the bootstraped parameters, and it will differ from
##' the RL computed with the estimated parameters in \code{object}.
##'
##' @param trace Integer level of verbosity.
##'
##' @param ... Arguments to be passed to \code{bs}.
##'
##' @return A data frame or an array with the RL and confidence limits
##' for each combination of \emph{date}, \emph{period} and
##' \emph{confidence level}. So when \code{confintMethod} is not
##' \code{"none"} the three values \code{Quant}, \code{L} and \code{U}
##' are given for each combination (the first not depending on the
##' confidence level). The data frame is in 'long' format, with
##'  different rows for multiple levels.
##'
##' @section Caution: When \code{confintMethod} is set to
##' \code{"loglik"} the required computing time is huge because
##' several constrained optimisations are required. Since the marginal
##' distribution hence the RL vary only slowly in practice, an
##' unnecessary computing burden will result when many values of
##' \code{newdate} are used as it is the case by default.
##'
##' @note Despite of its positive sounding name, the "bias correction"
##' can have a negative impact with some Extreme Value models as used
##' here, especially for non-parametric bootstap.
##' 
##' @author Yves Deville
##'
##' @examples
##' example(TVGEV)
##' p1 <- predict(res2, confintMethod = "none")
##' p2 <- predict(res2, confintMethod = "delta")
##' p3 <- predict(res2, confintMethod = "boot")
##' p4 <- predict(res2, newdate = "2020-01-01", confintMethod = "proflik") 
predict.TVGEV <- function(object,
                          newdate = NULL,
                          period = NULL,
                          level = 0.95,
                          confintMethod = c("delta", "none", "boot", "proflik"),
                          out = c("data.frame", "array"),
                          biasCorrect = FALSE,
                          trace = 1L,
                          ...) {

    if (!all(object$isCst)) {
        if (is.null(newdate)) {
            message("Since 'object' is really time-varying, the Return Levels\n",
                    "depend on the date. A default choice of dates is made here.\n",
                    "Use the 'newdate' formal to change this.")
        }
    } 
    
    out <- match.arg(out)
    parNames.GEV <- c("loc", "scale", "shape")
    
    if (is.null(period)) {
        period <- c(5, 10, 20, 50, 100, 150, 200, 500, 1000)
    } else {
        period <- sort(period)
    }

    nPeriod <- length(period)
    prob <- 1.0 - 1.0 / period
    fPeriod <- format(period)
    
    method <- match.arg(confintMethod)
    if (!missing(biasCorrect) && (method != "boot")) {
        warning("the argument 'biasCorrect' provided will not ",
                "be used since 'confintMethod' != \"boot\"")
    }
        
    ## take into account the order. 
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    if (is.null(newdate)) {
        newdate <- selectDate(object$data[ , object$date])
    }

    ## OLD CODE working with else            
    ## X <- object$X
    ## n <- object$n
    ## ind <- object$ind
    ## newdate <- object$data[ , object$date]
    ## fDate <- object$fDate

    n <- length(newdate)
    
    ## special case for non TV model
    if (!all(object$isCst)) {
        L <- modelMatrices.TSGEV(object, date = newdate)
        X <- L$X
    } else X <- NULL

    fDate <- format(newdate)
    
    theta <- psi2theta(object, date = newdate, deriv = TRUE)

    if (method == "none") {

        RL <- array(NA, dim = c(Date = n, Period = nPeriod),
                    dimnames = list(Date = fDate, Period = fPeriod))
        diagno <- NULL
        
        for (iPer in seq_along(period)) {
            RL[ , iPer] <- qGEV(p = prob[iPer],
                                loc = theta[ , 1],
                                scale = theta[ , 2],
                                shape = theta[ , 3])
        }
       
        
    } else if (method == "delta") {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        psiHat <- object$estimates
        covPsi <- vcov(object)
        p <- object$p
        
        if (is.null(covPsi) || any(is.na(covPsi))) {
            stop("vcov(object) must be a matrix with no NA. ",
                 "Consider changing 'confintMet'")
        }
        q <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)

        Quant <- array(NA,
                       dim = c(Date = n, Period = nPeriod),
                       dimnames = list(fDate, fPeriod))
        
        RL <- array(NA,
                    dim = c(Date = n, Period = nPeriod, Lim = 3L, Level = nLevel),
                    dimnames = list(Date = fDate, Period = fPeriod,
                        Type = c("Quant", "L", "U"), Level = fLevel)) 
        
        gradpsi <- array(NA, dim = c(n, p),
                         dimnames = list(fDate, object$parNames))

        diagno <- NULL
        
        for (iPer in seq_along(period)) {
 
            Quant[ , iPer] <- quant <- qGEV(p = prob[iPer],
                                         loc = theta[ , 1],
                                         scale = theta[ , 2],
                                         shape = theta[ , 3],
                                         deriv = TRUE)
            ## this is a n x p matrix
            gradtheta <- attr(quant, "gradient")
            
            for (nm in parNames.GEV) {
                if (!object$isCst[nm]) {
                    gradpsi[ , object$ind[[nm]]] <- 
                        sweep(X[[nm]], MARGIN = 1,
                              STATS = gradtheta[ , nm, drop = FALSE],
                              FUN = "*")
                } else {
                    gradpsi[ , object$ind[[nm]]] <- gradtheta[ , nm, drop = FALSE]
                }
            }

            ## standard deviation of the RL estimate.
            ## XXX could be faster with a 'crossprod' using the Cholesky root
            ## of 'covPsi'
            sdRL <- apply(gradpsi, 1, function(x)  t(x) %*% covPsi %*% x)

            ## could probably be done in one step?
            RL[ , iPer, "Quant",  ] <- Quant[ , iPer] 
            RL[ , iPer, "L",  ] <- Quant[ , iPer] + outer(sdRL, q[ , 1L]) 
            RL[ , iPer, "U",  ] <- Quant[ , iPer] + outer(sdRL, q[ , 2L])
            
        }
        
        ## change dim order: "Date", then "Level" then "Period" then "L or U"
        ## RL <- aperm(a = RL, perm = c(1, 4, 2, 3))
        
    } else if (method == "boot") {

        probL <- (1 - level) / 2
        probU <- 1 - probL
        
        if (is.null(object$boot)) {
            if (trace) {
                cat("'object' does not embed a bootstrap distribution.",
                    "Using `bs` to compute this.\n")
                
            }
            
            object$boot <- bs(object, ...)
            
        } else{
            message("'object' embeds a bootstrap distribution. Formal arguments\n",
                    " intented to be passed to `bs` (if any) will be ignored")
        }
        
        RL <- array(NA,
                    dim = c(Date = n, Period = nPeriod, Lim = 3L, Level = nLevel),
                    dimnames = list(Date = fDate, Period = fPeriod,
                        Type = c("Quant", "L", "U"), Level = fLevel)) 

        diagno <- NULL
        
        myFun <- function(psi, prob) {
            theta <- psi2theta(model = object, psi = psi, date = newdate)
            qGEV(prob, theta[ , 1], theta[ , 2], theta[ , 3])
        }
        
        for (iPer in seq_along(period)) {
            ## 'Quant' is a matrix with n rows and B columns. XXX Caution if n == 1.
            ## the dimension is lost!
            Quant <- apply(object$boot$estimate, MARGIN = 1, FUN = myFun, prob = prob[iPer])
            if (n == 1L) {
                dim(Quant) <- c(n, nrow(object$boot$estimate))
            }
            RL[ , iPer, "L", ] <- t(apply(Quant, MARGIN = 1, FUN = quantile, prob = probL))
            RL[ , iPer, "U", ] <- t(apply(Quant, MARGIN = 1, FUN = quantile, prob = probU))
            if (biasCorrect) {
                RL[ , iPer, "Quant", ] <- t(apply(Quant, MARGIN = 1, mean))
            } else {
                RL[ , iPer, "Quant", ] <- rep(qGEV(p = prob[iPer],
                                                   loc = theta[ , 1],
                                                   scale = theta[ , 2],
                                                   shape = theta[ , 3],
                                                   deriv = TRUE), times = nLevel)
            }
        }
        
    } else if (method == "proflik") {
        
        constrCheck <- -5e-3
        
        psiHat <- object$estimates
        
        RL <- array(NA,
                    dim = c(Date = n, Period = nPeriod, Lim = 3L, Level = nLevel),
                    dimnames = list(Date = fDate, Period = fPeriod,
                        Type = c("Quant", "L", "U"), Level = fLevel))

        diagno <- array(NA,
                        dim = c(Date = n, Period = nPeriod, Lim = 2L, Level = nLevel,
                                Diag = 3L),
                        dimnames = list(Date = fDate, Period = fPeriod,
                            Type = c("L", "U"), Level = fLevel,
                            Diag = c("status", "constraint", "gradDist")))
        
        ## ===================================================================
        ## For each parameter, we maximise/minimise it under the constraint
        ## that the logLik remains >= max logLik - delta where delta :=
        ## qchisq(1 - alpha) where alpha is given by the cofidence level.
        ##
        ## ===================================================================

        opts1 <- list("algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-5,
                      "ftol_abs" = 1.0e-7, "ftol_rel" = 1.0e-5,
                      "maxeval" = 3000,
                      "check_derivatives" = FALSE,
                      ## "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-8,
                      ##    "maxeval" = 3000,
                      ##    "ftol_abs" = 1.0e-4,
                      ##    "ftol_rel" = 1.0e-8),
                      "print_level" = 0)

        if (trace >= 2) {
            opts1[["check_derivatives"]] <- TRUE
            opts1[["check_derivatives_print"]] <- "all"
        }
         
         ## ==============================================================
         ## note that some arguments such as 'level' are unused but are
         ## required by the constraint
         ## ===============================================================
         
         f <- function(psi, level, prob, iDate, chgSign = FALSE, object) {

             theta <- psi2theta(model = object,
                                psi = psi,
                                date = newdate[iDate],
                                deriv = TRUE, checkNames = FALSE)
             
             RL <- qGEV(prob, theta[ , 1], theta[ , 2], theta[ , 3],
                        deriv = TRUE)

             ## 'nloptr' fails on NA and NaN!
             if (is.na(RL)) {
                 if (chgSign) {
                     return(list("objective" = Inf,
                                 "gradient" = rep(NaN, object$p)))
                 } else {
                     return(list("objective" = Inf,
                                 "gradient" = rep(NaN, object$p)))
                 }
             }
             
             gradtheta <- attr(RL, "gradient")

             if (all(object$isCst)) {
                 gradpsi <- gradtheta
                 names(gradpsi) <- object$parNames
             } else {
                 gradpsi <- rep(1.0, object$p)
                 names(gradpsi) <- object$parNames
                 
                 for (nm in parNames.GEV) {
                     if (!object$isCst[nm]) {
                         gradpsi[object$ind[[nm]]] <-
                             gradtheta[ , nm] * X[[nm]][iDate, ]
                     } else {
                         gradpsi[object$ind[[nm]]] <- gradtheta[1, nm]
                     }
                 }
             }

             
             if (chgSign) {
                 return(list("objective" = -RL, "gradient" = -gradpsi))
             } else {
                 return(list("objective" = RL, "gradient" = gradpsi))
             }
         }
         
         ## ==============================================================
         ## Up to (possibly unused) arguments, the constraint function
         ## is the same as that used in the 'confint' method.
         ## ===============================================================
         
         g <- function(psi, level, prob, iDate, chgSign = FALSE, object) {
             
             ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
             res <- object$negLogLikFun(psi = psi,
                                        deriv = TRUE,
                                        object = object)
             list("constraints" = res$objective - ellL,
                  "jacobian" = res$gradient)
             
         }
         
         ## =====================================================================
         ## note that although we recompute the gradient of the
         ## objective and the quantile of the chi-square distribution,
         ## this might be faster than re-defining the functions in the
         ## loop. Some experimentations would be needed to confirm
         ## this.
         ## =====================================================================

         for (iDate in 1L:n) { 

             ## if (trace) cat(sprintf("\no Finding CI for date %s\n", fDate[iDate]))
                
             for (iPer in seq_along(period)) {
                 
                 quant <- qGEV(p = prob[iPer],
                               loc = theta[iDate, 1],
                               scale = theta[iDate, 2],
                               shape = theta[iDate, 3],
                               deriv = FALSE)
                 
                ##  if (trace) cat(sprintf("\n   o period %d\n   ============\n",
                ##                         period[iPer]))
                 
                 for (iLev in rev(seq_along(level))) {
                     
                     lev <- level[iLev]

                     RL[iDate, iPer, "Quant", iLev] <- quant
                 }
             }
         }

        ## =========================================================
        ## Lower and upper bounds
        ## ==========================================================
        
        labs <- c("L" = "Lower", "U" = "Upper")
        sign <- c("L" = 1.0, "U" = -1.0)
        chgSign <- c("L" = 0.0, "U" = 1.0)
        
        if (trace) cat("\no Finding CI for Return Levels\n")

        for (LU in c("L", "U")) {
            
            if (trace) {
                cat(sprintf("\n**************\n %s bounds \n**************\n",
                            labs[LU]))
            }
            
            for (iLev in rev(seq_along(level))) {
                
                lev <- level[iLev]

                if (trace) {
                    cat(sprintf("Confidence Level: %5.2f\n", lev))
                }
                
                for (iDate in 1L:n) { 
                    
                    if (trace) cat(sprintf("\n  o Date: %s\n", fDate[iDate]))
                    
                    for (iPer in seq_along(period)) {

                        if (trace) cat(sprintf("\n    - Period:  %5d\n", period[iPer]))
                        
                        if ((iPer > 1L) && !is.null(psiIniPrec)) {
                            psi0 <- psiIniPrec
                        } else {
                            psi0 <- psiHat
                        }
                        
                        resOpt <- try(nloptr::nloptr(x0 = psi0,
                                                     eval_f = f,
                                                     eval_g_ineq = g,
                                                     level = lev,
                                                     prob = prob[iPer],
                                                     iDate = as.double(iDate),
                                                     chgSign = chgSign[LU],
                                                     opts = opts1,
                                                     object = object))
                        
                        diagno[iDate, iPer, LU, iLev, "status"] <-
                            resOpt$status
                        
                        if (trace == 1L) {
                            cat(sprintf("        Optimisation status: %d\n", resOpt[["status"]]))
                            names(resOpt$solution) <- object$parNames
                            cat(sprintf("        Objective: %7.2f\n", resOpt[["objective"]]))
                        } else  if (trace > 2L) {
                            cat("\nSOLUTION\n")
                            print(resOpt)
                        }
                        
                        ## The constraint must be active. We have to check that!
                        checkg <- g(psi = resOpt$solution,
                                    level = lev,
                                    prob = prob[iPer],
                                    iDate = iDate,
                                    chgSign = chgSign[LU],
                                    object = object)
                        
                        diagno[iDate, iPer, LU, iLev, "constraint"] <-
                            checkg$constraints
                        
                        ## The gradient of the objective must be colinear to the
                        ## jacobian of the constraint. We have to check that!
                        checkf <- f(psi = resOpt$solution,
                                    level = lev,
                                    prob = prob[iPer],
                                    iDate = iDate,
                                    chgSign = chgSign[LU],
                                    object = object)
                        
                        gradDist <- distLines(x1 = checkg$jacobian,
                                              x2 = checkf$gradient)
                        
                        diagno[iDate, iPer, LU, iLev, "gradDist"] <- gradDist
                        
                        if (trace) {
                            cat(sprintf("        Constraint check %10.7f\n",
                                        checkg$constraints))
                            cat(sprintf("        Gradient directions: %7.4f\n", gradDist))
                            if (trace == 2) {
                                cat("        Gradients\n")
                                print(rbind("        g" = checkg$jacobian,
                                            "        f" = checkf$gradient))
                                cat("\n")
                            }
                        }
                        check <- (checkg$constraints > constrCheck)
                        
                        if (!inherits(resOpt, "try-error") && (resOpt$status >= 0) && check) {
                            psiIniPrec <- resOpt[["solution"]]
                            RL[iDate, iPer, LU, iLev] <- sign[LU] * resOpt[["objective"]]
                        } else {
                            psiIniPrec <- NULL
                        }
                    }
                }
            }
        }
  

    }
    
    if (out == "data.frame") {
        
        if (requireNamespace("reshape2", quietly = TRUE)) {
            
            if (method == "none") {
                RL <- melt(RL, value.name = "Quant", varnames = c("Date", "Period"))
                RL$Date <- as.Date(RL$Date)
            } else {
                ## ====================================================================
                ## UGGLY CODE: there must be a simpler and more
                ## efficient way of doing that. The problem is that we
                ## want to drop the " Type" dimension but not the
                ## "Level" dimension even when it has no extension
                ## ====================================================================
                
                df <- list()
                for (nm in c("Quant", "L", "U")) {
                    RL1 <- RL[ , , nm, , drop = FALSE]
                    df[[nm]]  <- melt(RL1,
                                      value.name = nm,
                                      varnames = c("Date", "Period", "Type", "Level"))
                }
                RL <- data.frame(df[["Quant"]][, c("Date", "Period", "Level", "Quant")],
                                 L = df[["L"]]$L, U = df[["U"]]$U)
                
                RL$Date <- as.Date(RL$Date)
                ind <- with(RL, order(Date, Level, Period))
                RL <- RL[ind, , drop = FALSE]
            }
            
        } else {
            stop("the package 'reshape2' could not be used")
        }

        class(RL) <- c("predict.TVGEV", "data.frame")
    } 

    attr(RL, "diagno") <- diagno
    
    return(RL)

    
}


## ***********************************************************************
##' Plot predict Results for \code{TVGEV}.
##' 
##' @title Plot Predict Results for \code{TVGEV}
##'
##' @param x An object with class \code{"predict.TVGEV"} as returned
##' by the \code{predict} method. This must be a data frame so the
##' \code{predict} method must be called without using \code{out} or
##' setting it to the default value \code{"data.frame"}.
##' 
##' @param y Not used.
##'
##' @param gg Logical. If \code{TRUE} the plot is produced with the
##' \bold{ggplot2} package.
##' 
##' @param ... Not used yet.
##'
##' @return An object with class \code{"gg"}. Can be used with the method
##' \code{plot}, or equivalently with \code{print}.
##'
##' @note This function is intended to work with predictions computed
##' for a possibly large number of periods (to get smooth curves) but
##' with only a small number of dates, each apearing in a facet. So
##' the number of dates is limited to \code{6}.
##' 
##' @examples
##' example(TVGEV)
##' pred <- predict(res2, newdate = c("1960-01-01", "2000-01-01", "2020-01-01"),
##'                 level = c(0.70, 0.95), confintMethod = "delta")
##' g <- plot(pred)
plot.predict.TVGEV <- function(x, y, gg = TRUE, ... ) {

    ## avoid NOTE 'no visible binding' in checks
    Period <- L <- U <- Level <- Quant <- NULL    
    
    if (!gg) stop("Sorry, but for now only the 'gg = TRUE' is possible!")
    
    if (!requireNamespace("ggplot2", quietly = TRUE))  stop("package 'ggplot2' needed here")
        
    confLev <- attr(x, "confLevel")
    
    if (nd <- length(unique(x$Date)) > 6L) {
        stop(nd, "dates found in 'x'. The maximum allowed is 6")
    } 
    
    fill <- "darkgray"
    ## Find out if the confidence levels
    g1 <- ggplot(data = x)

    if (!is.null(x$L) && !is.null(x$U)) {
        g1 <- g1  + geom_ribbon(mapping = aes(x = Period, ymin = L, ymax = U, group = Level,
                                    ## colour = Level,
                                    fill = Level),
                                ## group = type,
                                alpha = 0.2)
        
        g1 <- g1 + geom_line(mapping = aes(x = Period, y = L, group = Level), alpha = 0.2)
        g1 <- g1 + geom_line(mapping = aes(x = Period, y = U, group = Level), alpha = 0.2)
    }
    ## g1 <- g1 + scale_fill_manual(name = "legend_name", values = c("red", "green"),
    ##                               labels = c('one','two'))
    g1 <- g1 + geom_line(data = x,
                         mapping = aes(x = Period, y = Quant, colour = "orangered"), size = 1)
    g1 <- g1 + theme_bw()
    g1 <- g1 + scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000))
    ## g
    g1 <- g1 + theme(axis.text.x = element_text(angle = 0, hjust = 0),
                     axis.title = element_text(face = 'bold', size = 12),
                     title = element_text(face = 'bold', size = 12),
                     legend.position = 'right',
                     legend.title = element_blank(),
                     legend.text = element_text(size = 12))

    g1 <- g1 + scale_colour_manual(name = "", values = "orangered", labels = "Quantile")

    if (!is.null(x$L) && !is.null(x$U)) {
        g1 <- g1 + scale_fill_manual(name = "Level", values = c("darkgray", "darkgray"))
    }
    
    g1 <- g1 + xlab("Period") + ylab("Quantile") 
    
    g1 <- g1 + facet_wrap( ~ Date)
        
    g1
   
}

