##*****************************************************************************
##' Generalised Residuals for a TVGEV model.
##' 
##' @title Generalised Residuals for a TVGEV Model
##'
##' @aliases resid.TVGEV
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param type The approximate distribution wanted.
##'
##' @param ... Not used yet.
##' 
##' @return A vector of generalised residuals which should
##' \emph{approximately} be independent and \emph{approximately}
##' follow the standard exponential or the uniform distribution,
##' depending on the value of \code{type}.
##' 
##' @note The upper 95\% quantile of the standard exponential is close
##' to \eqn{3} which can be used to gauge "large residuals".
##'
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' tv <- TVGEV(data = df, response = "TXMax", date = "Date",
##'             design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'             loc = ~ t1 + t1_1970)
##' e <- resid(tv)
##' plot(e)
##' mu <- tv$theta[ , "loc"]
##' plot(mu, e, type = "p", pch = 16, col = "darkcyan",
##'      main = "generalised residuals against 'loc'")
residuals.TVGEV <- function(object,
                            type = c("exp", "unif"),
                            ...) {

    type <- match.arg(type)
    Y <- object$data[ , object$response]
    theta <- psi2theta(model = object, psi = NULL)
  
    if (type == "exp") {
        e <- pGEV(Y, loc = theta[, 1L], scale = theta[, 2L], 
                  shape = theta[, 3L], lower.tail = FALSE)
        e <- -log(e)
    } else {
        e <- pGEV(Y, loc = theta[, 1L], scale = theta[, 2L], 
              shape = theta[, 3L], lower.tail = TRUE)
    }
    names(e) <- rownames(theta)
    attr(e, "date") <- object$data[ , object$date]
    attr(e, "type") <- type
    class(e) <- "resid.TVGEV"
    e
    
}

##*****************************************************************************
##' Plot the residuals of a \code{TVGEV} model object against the
##' date.
##'
##' @title Plot Residuals of a \code{TGVEV} model.
##'
##' @param x An object with class \code{"TVGEV"}.
##'
##' @param y Not used yet.
##'
##' @param ... Further arguments to be passed to \code{plot}.
##'
##' @return Nothing.
##'
##' @seealso \code{\link{residuals.TVGEV}}.
##' 
plot.resid.TVGEV <- function(x, y = NULL, ...) {

    if (!missing(y) && !is.null(y)) {
        warning("'y' formal not used by this method")
    }
    type <- attr(x, "type")
    plot(attr(x, "date"), x, type = "o",
         pch = 16, col = "orangered",
         xlab = "date",
         ylab = sprintf("residual, type = \"%s\"", type),
         ...)

    lims <- c(0.025, 0.975)

    if (type == "exp") lims <- -log(1 - lims)
    abline(h = lims, col = "SpringGreen3")
    
}

##*****************************************************************************
##' Confidence intervals for a \code{TVGEV} object
##'
##' @title Confidence Intervals for a \code{TVGEV} Object
##' 
##' @param object An object with class \code{"TVGEV"}.
##'
##' @param parm Parameter name. NOT USED YET.
##' 
##' @param level Confidence level.
##' 
##' @param method \code{"delta"}, \code{"proflik"},
##' \code{"boot"}.
##' 
##' @param trace Integer level of verbosity.
##' 
##' @param round Logical. If \code{TRUE} the confidence limits will be
##' rounded to a small number of digits. This number is chosen using the
##' smallest of the standard deviations for the estimated parameters.
##'
##' @param out Character giving the class of the output. By default a
##' three-dimensional array with dimensions: \emph{parameter},
##' \emph{lower/upper} limit, and \emph{level}. If \code{level} has length
##' \code{1}, using \code{drop} on the output will remove the third dimension
##' and produce a matrix. The \code{"data.frame"} gives the same results
##' in \code{'long'} format.
##' 
##' @param ... Arguments to be passed to the \code{profLik} or to the
##' \code{bs} method. For instance, they can be used to chose the type
##' of bootstrap or the number of bootstrap replications if
##' \code{method} is \code{"boot"}.
##'
##' @return An array or a data frame with the lower and upper bounds
##' \code{"L"} and \code{"U"} of the confidence intervals.
##'
##' @note For the bootstrap method(s), the time required depends on
##' whether \code{object} embeds a bootstrap distribution or not. When
##' no bootstrap distribution is found, it has to be computed using
##' the \code{bs} method; the formal arguments that are given in
##' \code{\dots} and intended to be used by \code{bs} will then be
##' ignored with a warning.
##'
##' @examples
##'
##' ## use the example of TVGEV which defines a TVGEV object named
##' ## 'res2'
##' example(TVGEV)
##' 
##' ## (default) delta method, several confidence levels
##' ci <- confint(res2, level = c(0.90, 0.70))
##' ci1 <- confint(res2, level = c(0.95, 0.90, 0.70), out = "data.frame")
##' ci2 <- confint(res2, level = 0.95, out = "data.frame")
##'
##' ## Profile-likelihood
##' ci3 <- confint(res2, level = 0.95, method = "proflik")
##' 
confint.TVGEV <- function(object,
                          parm = NULL, 
                          level = 0.95,
                          method = c("delta", "boot", "proflik"),
                          trace = 1L,
                          round = TRUE,
                          out = c("array", "data.frame"),
                          ...) {

    out <- match.arg(out)
    
    method <- match.arg(method)

    ## take into account the order. 
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    
    if (method == "delta") {
        
        probL <- (1 - level) / 2
        probU <- 1 - probL
        psiHat <- object$estimates
        sigHat <- object$sd
        q <- qnorm(cbind(probL, probU), mean = 0.0, sd = 1.0)

        ci <- array(psiHat, dim = c(object$p, 2L, nLevel),
                   dimnames = list(object$parNames, c("L", "U"), fLevel)) 
        cw <-  array(sigHat, dim = c(object$p, 2L, nLevel),
                     dimnames = list(object$parNames, c("L", "U"), fLevel)) 

        cw <- sweep(cw, MARGIN = c(3, 2), STATS = q, FUN = "*")
        ci <- ci + cw
        
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
        
        ci <- apply(object$boot$estimates, 2, quantile,
                    prob = rbind(probL, probU))
        ci <- array(ci, dim = c(2L, nLevel, object$p),
                     dimnames = list(c("L", "U"), fLevel, object$parNames))
        ci <- aperm(ci, perm = c(3L, 1L, 2L))
        
    } else if (method == "proflik") {

        prob <- 1 - level
        psiHat <- object$estimates
        ci <- array(NA, dim = c(object$p, 2L, nLevel),
                    dimnames = list(object$parNames, c("L", "U"), fLevel)) 

        ## ===================================================================
        ## For each parameter, we maximise/minimise it under the constraint
        ## that the logLik remains >= max logLik - delta where delta :=
        ## qchisq(1 - alpha) where alpha is given by the cofidence level.
        ##
        ## ===================================================================

        opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel" = 1.0e-10,
                      "maxeval" = 1000,
                      ## "check_derivatives" = TRUE, "check_derivatives_print" = "all",
                      "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7),
                      "print_level" = 0)

        ## note that some arguments sucha as 'level' is unused but is
        ## required by the constraint
        
        f <- function(psi, k, chgSign = FALSE, level, object) {
            grad <- rep(0.0, object$p)
            if (chgSign) {
                grad[k] <- -1.0
                return(list("objective" = -psi[k], "gradient" = grad))
            } else {
                grad[k] <- 1.0
                return(list("objective" = psi[k], "gradient" = grad))
            }
        }

        g <- function(psi, k, chgSign = FALSE, level, object) {
            
            ellL  <-  object$negLogLik + qchisq(level, df = 1) / 2.0
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
       
        for (k in 1L:object$p) {
            
            if (trace) cat(sprintf("o Finding CI for \"%s\n", object$parNames[k]))

            ilevPrec <- 1L
            
            for (ilev in rev(seq_along(level))) {
                
                lev <- level[ilev]
                
                if (trace) {
                    cat(sprintf("  %s, lower bound\n", fLevel[ilev]))
                }

                ## =========================================================
                ## if we have successfully computed the result for a
                ## larger confidence level (and the same parameter),
                ## use it as initial guess
                ## ==========================================================
                
                if ((ilevPrec > 1L) && (!is.null(psiLPrec))) {
                    psi0 <- psiLPrec
                } else {
                    psi0 <- psiHat
                }
                
                resL <- try(nloptr::nloptr(x0 = psi0,
                                           eval_f = f,
                                           eval_g_ineq = g,
                                           k = k, level = lev, chgSign = FALSE,
                                           opts = opts1,
                                           object = object))
                            
                if (trace > 1L) {
                    cat("SOLUTION\n")
                    print(resL)
                }

                
                ## the constraint must be active
                check <- g(resL$solution, object = object, k = k, level = lev)$constraints
                check <- (check < 1e-5)
                
                if (!inherits(resL, "try-error") && (resL$status >= 0) && check) {
                    psiLPrec <- resL[["solution"]]
                    ci[k, "L", ilev] <- psiLPrec[k]
                } else {
                    psiLPrec <- NULL
                }
                
                ## here we maximise will 'nloptr' only minimises things
                if (trace) {
                    cat(sprintf("  %s, upper bound\n", fLevel[ilev]))
                }

                ## =========================================================
                ## if we have successfully computed the result for a
                ## larger confidence level (and the same parameter),
                ## use it as initial guess
                ## ==========================================================
                
                if ((ilevPrec > 1L) && (!is.null(psiUPrec))) {
                    psi0 <- psiUPrec
                } else {
                    psi0 <- psiHat
                }
            
                resU <- try(nloptr::nloptr(x0 = psi0,
                                           eval_f = f,
                                           eval_g_ineq = g,
                                           k = k, level = lev, chgSign = TRUE,
                                           opts = opts1,
                                           object = object))
                if (trace > 1L) {
                    cat("SOLUTION\n")
                    print(resU)
                }
                
                ## the constraint must be active
                check <- g(resU$solution, object = object, k = k, level = lev)$constraints
                check <- (check < 1e-5)
             
                if (!inherits(resU, "try-error") && (resU$status >= 0) && check) {
                    psiUPrec <- resU[["solution"]]
                    ci[k, "U", ilev] <- psiUPrec[k]
                } else {
                    psiUPrec <- NULL
                }
                
                ilevPrec <- ilev
            }
        }
   
    }

    if (round && (!is.null(object$sd)) && (!any(is.na(object$sd)))) {
        d <- ceiling(-log(min(object$sd), 10)) + 1
        ci <- round(ci, digits = d)
    }
    
    if (out == "data.frame") {
    
        ci <- array(ci, dim = c(object$p * nLevel, 2L),
                    dimnames = list(rep(object$parNames, nLevel), c("L", "U")))
        ci <- data.frame(parm = rep(object$parNames, times = nLevel),
                         level = rep(fLevel, each = object$p),
                         L = ci[ , "L"], U = ci[ , "U"])
    }
    
    ci 
}
