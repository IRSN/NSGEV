
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
##' @importFrom stats confint
##' @method confint TVGEV
##' @export
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
        psiHat <- object$estimate
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
        
        ci <- apply(object$boot$estimate, 2, quantile,
                    prob = rbind(probL, probU))
        ci <- array(ci, dim = c(2L, nLevel, object$p),
                     dimnames = list(c("L", "U"), fLevel, object$parNames))
        ci <- aperm(ci, perm = c(3L, 1L, 2L))
        
    } else if (method == "proflik") {

        prob <- 1 - level
        psiHat <- object$estimate
        ci <- array(NA, dim = c(object$p, 2L, nLevel),
                    dimnames = list(object$parNames, c("L", "U"), fLevel)) 

        ## ===================================================================
        ## For each parameter, we maximise/minimise it under the constraint
        ## that the logLik remains >= max logLik - delta where delta :=
        ## qchisq(1 - alpha) where alpha is given by the cofidence level.
        ##
        ## ===================================================================

        opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel" = 1.0e-8, "ftol_abs" = 1.0e-8, "ftol_rel" = 1.0e-8,
                      "maxeval" = 3000,
                      "check_derivatives" = FALSE,
                      "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-8,
                          "maxeval" = 3000,
                          "ftol_abs" = 1.0e-8,
                          "ftol_rel" = 1.0e-8),
                      "print_level" = 0)

        if (trace >= 2) {
            opts1[["check_derivatives"]] <- TRUE
            opts1[["check_derivatives_print"]] <-  "all"
        }
        if (trace >= 3) {
              opts1[["print_level"]] <- 1
        }

        ## ==============================================================
        ## note that some arguments such as 'level' are unused but are
        ## required by the constraint
        ## ==============================================================
        
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

            psi1 <- psi
            names(psi1) <- object$parNames
            sigma <- min(psi2theta(model = object, psi = psi1)[ , 2])

            ellL  <-  object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(psi = psi, deriv = TRUE, object = object)

            if (!is.na(sigma) && (sigma > 0)) {
                res2 <- list("constraints" = res$objective - ellL,
                             "jacobian" = res$gradient)
            } else {
                res2 <-  list("constraints" = 1,
                              "jacobian" = rep(NaN, object$p))
            }
            res2 
            
        }

        ## =====================================================================
        ## note that although we recompute the gradient of the
        ## objective and the quantile of the chi-square distribution,
        ## this might be faster than re-defining the functions in the
        ## loop. Some experimentations would be needed to confirm
        ## this.
        ## =====================================================================
       
        for (k in 1L:object$p) {
            
            if (trace) cat(sprintf("\n\no Finding CI for \"%s\"\n", object$parNames[k]))

            ilevPrec <- 1L
            
            for (ilev in seq_along(level)) {
                
                lev <- level[ilev]
                
                if (trace) {
                    cat(sprintf("\n   o %s, lower bound: ", fLevel[ilev]))
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
                if (trace == 1L) {
                    cat(sprintf("%7.2f\n", resL[["objective"]])) 
                } else if (trace > 1L) {
                    cat("SOLUTION\n")
                    print(resL)
                }

                
                ## the constraint must be active
                check <- g(resL$solution, object = object, k = k, level = lev)$constraints

                check2 <- object$negLogLikFun(psi = resL$solution,
                                              deriv = FALSE,
                                              object = object)

                if (trace) cat(sprintf("   Constraint check %10.7f, %10.4f\n", check, check2))

                check <- (check > -1e-3)
                
                if (!inherits(resL, "try-error") && (resL$status %in% c(3, 4))) {
                    psiLPrec <- resL[["solution"]]
                    ci[k, "L", ilev] <- psiLPrec[k]
                } else {
                    psiLPrec <- NULL
                }
                
                ## here we maximise will 'nloptr' only minimises things
                if (trace) {
                    cat(sprintf("   %s, upper bound: ", fLevel[ilev]))
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
                
                if (trace == 1L) {
                    cat(sprintf("%7.2f\n", -resU[["objective"]])) 
                } else  if (trace > 1L) {
                    cat("SOLUTION\n")
                    print(resU)
                }
                
                ## the constraint must be active
                check <- g(resU$solution, object = object, k = k, level = lev)$constraints

                check2 <- object$negLogLikFun(psi = resU$solution,deriv = FALSE, object = object)

                if (trace) cat(sprintf("   Constraint & value %10.7f, %10.4f\n", check, check2))

                check <- (check > -1e-3)
                
                if (!inherits(resU, "try-error") && (resU$status %in% c(3, 4))) {
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
