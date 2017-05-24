## *************************************************************************
##' Profile-likelihood inference method.
##'
##' @title Profile-Likelihood Inference Method
##'
##' @param object An object representing a parametric model.
##'
##' @param fun A numeric function of the vector of parameters.
##' 
##' @param ... Further arguments for methods.
##' 
##' @return The result, typically a numeric array.
##'
##' @author Yves Deville
##' 
profLik <- function(object, fun, ...) {
    UseMethod("profLik")
}

## *************************************************************************
##' Profile-likelihood inference for \code{TVGEV} objects.
##'
##' Compute the lower and upper end-points of a profile-likelihood
##' based confidence interval. The (apparently new) method used here
##' relies on maximising the function of interest, say
##' \eqn{\rho(\boldsymbol(\psi)}{\rho(\psi)}, under the constraint
##' that the log-likelihood is greater than the maximal log-likelihood
##' minus a positive quantity \eqn{\delta} depending on the confidence
##' level. This differs from the usual method which relies on a
##' univariate zero-finding for the profile-likelihood function, the
##' evaluation of which relies on a \eqn{p-1} dimensional
##' optimisation. As a major advantage, the new method does not require
##' a re-parameterisation of the model.
##' 
##' @title  Profile-Likelihood Inference for \code{TVGEV} Objects
##' 
##' @param object A \code{TSGEV} object.
##'
##' @param fun A function of the parameter for which the
##' profile-likelihood will be carried over. This function must have
##' the arguments: \code{psi} for the vector of parameters and
##' \code{object} for the model object. If needed, a wrapper function
##' can be used use more arguments, see \bold{Details}.
##' 
##' @param level Level of confidence. Can be of length \code{> 1}.
##'
##' @param deriv Logical. If \code{TRUE}, the function \code{fun} is
##' assumed to provide a gradient vector as an attribute named
##' \code{"gradient"} of the result. For now \code{deriv} can only be
##' \code{TRUE}.
##'
##' @param trace Level of verbosity; \code{trace = 0} prints nothing.
##'
##' @param ... Not used yet. 
##'
##' @return An array with the value of the function and the
##' corresponding Lower and Upper end-points for the given confidence
##' levels.
##'
##' @author Yves Deville
##'
##' @references
##'
##' Deville Y. (2017) "Profile-likelihood using constrained
##' optimisation". Unpublished Tech. Report.
##' 
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' ## fit a TVGEV model with constant parameters.
##' res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'               estim = "nloptr")
##' ## define a function of the parameter vector
##' myfun <- function(psi, object, deriv = TRUE) {
##'     res <- psi[1]
##'     grad <- rep(0.0, object$p)
##'     grad[1] <- 1
##'     attr(res, "gradient") <- grad
##'     res
##' }
##'
##' profLik.TVGEV(object = res1, fun = myfun, deriv = TRUE)
##' 
profLik.TVGEV <- function(object,
                          fun,
                          level = 0.95,
                          deriv = TRUE,
                          trace = 0,
                          ...) {

    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    psiHat <- object$estimates
    constrCheck <- -5e-3
    
    res <- array(NA, dim = c(Lim = 3L, Level = nLevel),
                 dimnames = list(Type = c("Quant", "L", "U"), Level = fLevel))

    ## ===================================================================
    ## For each parameter, we maximise/minimise it under the constraint
    ## that the logLik remains >= max logLik - delta where delta :=
    ## qchisq(1 - alpha) where alpha is given by the cofidence level.
    ##
    ## ===================================================================

    opts1 <- list("algorithm" = "NLOPT_LD_AUGLAG",
                  "xtol_rel" = 1.0e-8, "ftol_abs" = 1.0e-4,
                      "maxeval" = 3000,
                  "check_derivatives" = FALSE,
                  "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-8,
                      "maxeval" = 3000,
                      "ftol_abs" = 1.0e-4,
                      "ftol_rel" = 1.0e-8),
                  "print_level" = 0)
    
    if (trace >= 2) {
        opts1[["check_derivatives"]] <- TRUE
        opts1[["check_derivatives_print"]] <- "all"
    }
    
    if (deriv) {
        
        f <- function(psi, object, level, chgSign) {
            res <- fun(psi, object)
            gradpsi <-  attr(res, "gradient")
            if (chgSign) {
                 return(list("objective" = -res, "gradient" = -gradpsi))
             } else {
                 return(list("objective" = res, "gradient" = gradpsi))
             }
        }
    
        g <- function(psi, object, level, chgSign) {
            
            ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
            res <- object$negLogLikFun(psi = psi,
                                       object = object,
                                       deriv = TRUE,
                                       ...)
            list("constraints" = res$objective - ellL,
                 "jacobian" = res$gradient)
            
        }
        
    } else {
        
        stop("for now, the method is only implemented for deriv = TRUE")

    }
    
    for (iLev in rev(seq_along(level))) {
        lev <- level[iLev]
         if (trace) {
             cat(sprintf("     %s, lower bound: ", fLevel[iLev]))
         }
                     
        ## =========================================================
        ## if we have successfully computed the result for a
        ## larger confidence level (and the same parameter),
        ## use it as initial guess
        ## ==========================================================
        
        if ((iLev > 1L) && !is.null(psiLPrec)) {
            psi0 <- psiLPrec
        } else {
            psi0 <- psiHat
        }
        
        resL <- try(nloptr::nloptr(x0 = psi0,
                                   eval_f = f,
                                   eval_g_ineq = g,
                                   level = lev,
                                   chgSign = as.double(FALSE),
                                   opts = opts1,
                                   object = object))
        
        if (trace == 1L) {
            names(resL$solution) <- object$parNames
            cat(sprintf("%7.2f\n", resL[["objective"]]))
        } else  if (trace > 1L) {
            cat("\nSOLUTION\n")
            print(resL)
        }
        
        ## the constraint must be active: check that!
        check <- g(psi = resL$solution,
                   level = lev,
                   chgSign = FALSE,
                   object = object)$constraints
        
        check2 <- object$negLogLikFun(psi = resL$solution,
                                      deriv = FALSE,
                                      object = object)
        
        if (trace) cat(sprintf("     Constraint check %10.7f, %10.4f\n", check, check2))
        
        check <- (check > constrCheck)
        
        if (!inherits(resL, "try-error") && (resL$status >= 0) && check) {
            psiLPrec <- resL[["solution"]]
            res["L", iLev] <- resL[["objective"]]
        } else {
            psiLPrec <- NULL
        }
        
        ## here we maximise will 'nloptr' only minimises things
        if (trace) {
            cat(sprintf("     %s, upper bound: ", fLevel[iLev]))
        }
        
        ## =========================================================
        ## if we have successfully computed the result for a
        ## larger confidence level (and the same parameter),
        ## use it as initial guess
        ## ==========================================================
        
        if ((iLev > 1L) && !is.null(psiUPrec)) {
            psi0 <- psiUPrec
        } else {
            psi0 <- psiHat
        }
        
        resU <- try(nloptr::nloptr(x0 = psi0,
                                   eval_f = f,
                                   eval_g_ineq = g,
                                   level = lev,
                                   chgSign = as.double(TRUE),
                                   opts = opts1,
                                   object = object))
        
        if (trace == 1L) {
            names(resU$solution) <- object$parNames
            cat(sprintf("%7.2f\n", -resU[["objective"]]))
        } else if (trace > 1L) {
            cat("\nSOLUTION\n")
            print(resU)
        }
        
        ## the constraint must be active
        names(resU$solution) <- object$parNames
        check <- g(psi = resU$solution,
                   level = lev,
                   chgSign = TRUE,
                   object = object)$constraints
        
        check2 <- object$negLogLikFun(psi = resU$solution,
                                      object = object,
                                      deriv = FALSE)
                     
        if (trace) cat(sprintf("     Constraint check %10.7f, %10.4f\n", check, check2))
        check <- (check > constrCheck)
        
        if (!inherits(resU, "try-error") && (resU$status >= 0) && check) {
            psiUPrec <- resU[["solution"]]
            res["U", iLev] <- -resU[["objective"]]
        } else {
            psiUPrec <- NULL
        }
        
    }
    
    res

}
