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
##' @export
##' 
profLik <- function(object, fun, ...) {
    UseMethod("profLik")
}

## *************************************************************************
##' Profile-likelihood inference for \code{TVGEV} objects.
##'
##' Compute the lower and upper end-points of a profile-likelihood
##' based confidence interval. The (apparently new) method used here
##' relies on maximising and minimising the function of interest, say
##' \eqn{\eta(\boldsymbol(\psi)}{\eta(\psi)}, under the constraint
##' that the log-likelihood is greater than the maximal log-likelihood
##' minus a positive quantity \eqn{\delta} depending on the confidence
##' level. This differs from the usual method which relies on an
##' univariate zero-finding for the profile-likelihood function (minus
##' a constant). Remind that each evaluation of the profile requires a
##' \eqn{p-1} dimensional optimisation. As a major advantage, the new
##' method does not require a re-parameterisation of the model.
##' 
##' @title  Profile-Likelihood Inference for \code{TVGEV} Objects
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param fun A function of the parameter vector for which the
##' profile-likelihood will be carried over. This function must have
##' the arguments: \code{psi} for the vector of parameters and
##' \code{object} for the model object; so the function can use the
##' some of slots of \code{object}. If needed, a wrapper function can
##' be used use more arguments, see \bold{Details}.
##' 
##' @param level Level of confidence. Can be of length \code{> 1}.
##'
##' @param deriv Logical. If \code{TRUE}, the function \code{fun} is
##' assumed to provide a gradient vector as an attribute named
##' \code{"gradient"} of the result. For now \code{deriv} can only be
##' \code{TRUE}, which implies that \code{fun} \emph{must} compute the
##' gradient.
##'
##' @param trace Level of verbosity; \code{trace = 0} prints nothing.
##'
##' @param ... Not used yet. 
##'
##' @return An array with the value of the function and the
##' corresponding Lower and Upper end-points for the given confidence
##' levels. This array has two attributes with names \code{"diagno"}
##' and \code{"psi"} which both are arrays. The attributes provide
##' information about the numerical optimisation and the values of the
##' vector of parameter that maximised or minimised the function
##' \code{fun}.
##'
##' @author Yves Deville
##'
##' @references
##'
##' Deville Y. (2017) "Profile-likelihood using constrained
##' optimisation". Unpublished Tech. Report.
##'
##' @note For each confidence limit the numerical optimisation may
##' fail, in which case the limit will be \code{NA}. Using \code{trace
##' = 1} can be useful to further check the optimisation. The
##' \code{Optimisation status} value should be \code{3} or \code{4}
##' for small changes on the parameter or on the objective. On the
##' other hand, a value of \code{5} indicates that the maximal number
##' of iterations was reached, which is considered here as a
##' failure. The \code{Constaint check} value should be small because
##' the constraint must be active at the optimum. The \code{gradDist}
##' is the distance between the two directions of the gradient vectors
##' (objective and constraint). It should be small as well because the
##' gradients must be colinear at the optimum (Lagrange conditions).
##'
##' @method profLik TVGEV
##' @export
##' 
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##'
##' ## fit a TVGEV model with constant parameters.
##' res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'               estim = "nloptr")
##' 
##' ## define a function of the parameter vector: here the first component
##' ## This is for illustration only since the the result can be obtained
##' ## using the 'confint' method with \code{method = "proflik"}, which
##' ## gives the confidence intervals for each of the parameters.
##' 
##' myfun <- function(psi, object) {
##'     res <- psi[1]
##'     grad <- rep(0.0, object$p)
##'     grad[1] <- 1
##'     attr(res, "gradient") <- grad
##'     res
##' }
##'
##' pl <- profLik(object = res1, fun = myfun, deriv = TRUE)
##' 
##' confint(res1, method = "proflik")
profLik.TVGEV <- function(object,
                          fun,
                          level = 0.95,
                          deriv = TRUE,
                          trace = 0,
                          ...) {

    if (FALSE) {
        dots <- match.call(expand.dots = FALSE)[["..."]]
        nmOk <- names(dots) %in% names(formals(fun))
        
        if (!all(nmOk)) {
            stop("all formals passed through the dots '...' must be ",
                 "formals of the function given in 'fun'")
        }
    }
    
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    psiHat <- object$estimate
    constrCheck <- -5e-3
    
    res <- array(NA, dim = c(Lim = 3L, Level = nLevel),
                 dimnames = list(Type = c("est", "L", "U"), Level = fLevel))

    ## ===================================================================
    ## For each parameter, we maximise/minimise it under the constraint
    ## that the logLik remains >= max logLik - delta where delta :=
    ## qchisq(1 - alpha) where alpha is given by the cofidence level.
    ##
    ## The constrained optim is performed using an augmented
    ## Lagrangian method which requires a local companion algorithm
    ## with its own settings. So a sublist is used to tune the local
    ## optimisation.
    ## ===================================================================

    opts1 <-
        list("algorithm" = "NLOPT_LD_AUGLAG",
             "xtol_rel" = 1.0e-6, "ftol_abs" = 1.0e-6, "ftol_rel" = 1.0e-6,
             "maxeval" = 3000,
             "check_derivatives" = FALSE,
             "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                 "xtol_rel" = 1.0e-6,
                 "maxeval" = 3000,
                 "ftol_abs" = 1.0e-6,
                 "ftol_rel" = 1.0e-6),
             "print_level" = 0)
    
    if (trace >= 2) {
        opts1[["check_derivatives"]] <- TRUE
        opts1[["check_derivatives_print"]] <- "all"
    }

    ## CHNAGE on 2024-11-20. EXPERIMENTAL
    opts <- list()
    opts[[1]] <- list("algorithm" = "NLOPT_LD_MMA",
                      "xtol_rel" = 1.0e-5,
                      "ftol_abs" = 1.0e-7, "ftol_rel" = 1.0e-5,
                      "maxeval" = 200,
                      "check_derivatives" = FALSE,
                      "print_level" = 0)
    
    opts[[2]] <- list("algorithm" = "NLOPT_LD_AUGLAG",
                      "xtol_rel" = 1.0e-5,
                      "ftol_abs" = 1.0e-7, "ftol_rel" = 1.0e-5,
                      "maxeval" = 500,
                      "check_derivatives" = FALSE,
                      "local_opts" = list("algorithm" = "NLOPT_LD_MMA",
                                          "xtol_rel" = 1.0e-5,
                                          "maxeval" = 1000,
                                          "ftol_abs" = 1.0e-7,
                                          "ftol_rel" = 1.0e-5),
                      "print_level" = 0)

    
    if (!deriv) {
        stop("for now, the method is only implemented for deriv = TRUE")
    }

    f <- function(psi, object, level, chgSign) {

        ## DEBUG : it could be possible to get sigma <= 0 during optim
        ## sigma <- min(psi2theta(model = object, psi = psi, checkNames = FALSE)[ , 2])
        ## cat("sigma = ", sigma, "\n")

        res <- fun(psi, object)
        
        ## 'nloptr' fails on NA and NaN!
        if (is.na(res)) {
            if (chgSign) {
                return(list("objective" = Inf,
                            "gradient" = rep(NaN, object$p)))
            } else {
                return(list("objective" = Inf,
                            "gradient" = rep(NaN, object$p)))
            }
        }
        
        gradpsi <-  attr(res, "gradient")
        
        if (chgSign) {
            return(list("objective" = -res, "gradient" = -gradpsi))
        } else {
            return(list("objective" = res, "gradient" = gradpsi))
        }
    }
    
    ## changed on 2019-06-04: add a constraint on the minimum of
    ## the GEV scales which must be positive. Without this
    ## constraint, the optimisation accepts negative values of
    ## the scale and diverges.
    
    g <- function(psi, object, level, chgSign) {
        
        psi1 <- psi
        names(psi1) <- object$parNames
        sigma <- min(psi2theta(model = object, psi = psi1)[ , 2])
        
        ellL <- object$negLogLik + qchisq(level, df = 1) / 2.0
        res <- object$negLogLikFun(psi = psi, object = object, deriv = TRUE)
        
        if (!is.na(sigma) && (sigma > 0)) {
            res2 <- list("constraints" = res$objective - ellL,
                         "jacobian" = res$gradient)
        } else {
            res2 <-  list("constraints" = 1,
                          "jacobian" = rep(NaN, object$p))
        }
        res2 
    }
    
    ## the confidence level is not used here so we take any value
    val <- f(psiHat, object, level = 0.95, chgSign = FALSE)$objective
    res["est",  ] <- val
    
    ## ========================================================================
    ## keep some information about optimisation: diagnostics and value
    ## of the parameter vector which lead to the min or max of the
    ## profiled function.
    ## ========================================================================
    
    diagno <-
        array(NA,
              dim = c(Lim = 2L,
                  Level = nLevel,
                  Diag = 4L),
              dimnames = list(Type = c("L", "U"),
                  Level = fLevel,
                  Diag = c("status", "objective", "constraint", "gradDist")))
    
    Psi <-
        array(NA,
              dim = c(Lim = 2L,
                  Level = nLevel,
                  psi = object$p),
              dimnames = list(Type = c("L", "U"),
                  Level = fLevel,
                  psi = object$parNames))
    
    labs <- c("L" = "Lower", "U" = "Upper")
    sign <- c("L" = 1.0, "U" = -1.0)
    chgSign <- c("L" = 0.0, "U" = 1.0)
    
    for (LU in c("U", "L")) {
        
        for (iLev in seq_along(level)) {
            
            lev <- level[iLev]
            
            if (trace) {
                cat(sprintf("\n\no %s bound for level %s\n",
                            labs[LU], fLevel[iLev]))
            }
            
            ## =========================================================
            ## if we have successfully computed the result for a
            ## larger confidence level (and the same parameter), use
            ## the corresponding parameter as initial guess
            ## ==========================================================
            
            if ((iLev > 1L) && !is.null(psiPrec)) {
                psi0 <- psiPrec
                if (trace > 1) {
                    cat("\nInitialising with the previous conf. level\n")
                }
            } else {
                psi0 <- psiHat
                if (trace > 1) {
                    cat("\nInitialising with the MLE\n")
                }
                ## if (LU == "U") {
                ##     if (trace > 1) {
                ##         cat("\nInitialising with the MLE\n")
                ##     }
                ## } else {
                ##     if (!any(is.na(Psi["U", iLev, ]))) {
                ##         psi0 <- 2.0 * psi0 - 1.0 * Psi["U", iLev, ]
                ##         if (trace > 1) {
                ##             cat("\nInitialising by symmetrising the solution for the upper bound\n")
                ##         }
                ##     } else {
                ##         cat("\nInitialising with the MLE because of failure for \"U\"\n")
                ##     }
                ## }
            }
            
            resOpt <- try(nloptr::nloptr(x0 = psi0,
                                         eval_f = f,
                                         eval_g_ineq = g,
                                         level = lev,
                                         chgSign = chgSign[LU],
                                         opts = opts[[1]],
                                         object = object))

            diagno[LU, iLev, "status"] <- resOpt$status
            if (trace == 1L) {
                cat(sprintf("    Optimisation status: %d\n", resOpt$status))
                cat(sprintf("    Iterations:          %d\n", resOpt$iterations))
            }
            
            if (trace > 1L) {
                cat("\nSOLUTION\n")
                print(resOpt)
            }
            
            ## ================================================================
            ## compute value of the constaint as well as the distance
            ## between the two directions gradient of objective 'f'
            ## and gradient of constraint 'g' to see if the constaint
            ## is active at solution
            ## ================================================================

            if (!inherits(resOpt, "try-error") && (resOpt$status %in% c(3, 4))) {
                
                checkg <- g(psi = resOpt$solution,
                            level = lev,
                            chgSign = chgSign[LU],
                            object = object)
                
                checkf <- object$negLogLikFun(psi = resOpt$solution,
                                              deriv = TRUE,
                                              object = object)
                
                diagno[LU, iLev, "objective"] <- checkf$objective
                diagno[LU, iLev, "constraint"] <- checkg$constraints
                
                if (trace == 1L) {
                    cat(sprintf("    Objective value:  %10.7f\n", checkf$objective))
                    cat(sprintf("    Constraint check: %10.7f\n", checkg$constraints))
                }
                
                gradDist <- distLines(x1 = checkg$jacobian,
                                      x2 = checkf$gradient)
                
                diagno[LU, iLev, "gradDist"] <- gradDist
                
                if (trace == 1L) {
                    cat(sprintf("    gradDist:        %10.7f\n", gradDist))
                    ## print(rbind("    f " = checkf$gradient,
                    ##             "    g " = checkg$jacobian))
                    
                }

                if ( (!is.na(gradDist)) && (gradDist < 0.05)) {
                    optDone <- TRUE
                    psiPrec <- resOpt[["solution"]]
                    res[LU, iLev] <- sign[LU] * resOpt[["objective"]]
                    Psi[LU, iLev, ] <- resOpt[["solution"]]
                } else {
                    psiPrec <- NULL
                }
                
            } else {
                psiPrec <- NULL
            }
      
        }
        
    }

    ## attach diagnostic and parameter values as attributes.
    attr(res, "diagno") <- diagno
    attr(res, "psi") <- Psi
    attr(res, "class") <- "profLik.TVGEV"
    invisible(res)

}

##' @method print profLik.TVGEV
##' @export
##' @keywords internal
##' 
print.profLik.TVGEV <- function(x, diagno = FALSE, ...) {

    if (!diagno) {
        for (nm in c("diagno", "psi", "class")) {
            attr(x, nm) <- NULL
        }
        print(x, ...)
    } else {
        diagnoDat <- attr(x, "diagno")
        for (diagnm in dimnames(diagnoDat)[["Diag"]]) {
            cat(sprintf("o %s\n", diagnm))
            print(diagnoDat[ , , diagnm])
        }
    }
    
}

##' @method summary profLik.TVGEV
##' @export
##' @keywords internal
##' 
summary.profLik.TVGEV <- function(object, ...) {
    
    profLiked <- object
    for (nm in c("diagno", "psi", "class")) {
        attr(profLiked, nm) <- NULL
    }
    
    diagnoDat <- attr(object, "diagno")
    for (diagnm in dimnames(diagnoDat)[["Diag"]]) {
        cat(sprintf("o %s\n", diagnm))
        print(diagnoDat[ , , diagnm])
    }

    L <- list(profLiked = profLiked,
              diagnoDat = diagnoDat)
    class(L) <- "summary.profLik.TVGEV"
    L
}

##' @method print summary.profLik.TVGEV
##' @export
##' @keywords internal
##' 
print.summary.profLik.TVGEV <- function(x, ...) {
    cat("o Profiled quantity\n")
    print(unclass(x$profLiked))
    cat("\n")
    cat("o Optimisation diagnostics\n")
    for (diagnm in dimnames(x$diagnoDat)[["Diag"]]) {
        cat(sprintf("   o %s\n", diagnm))
        print(x$diagnoDat[ , , diagnm])
    }
    
}
