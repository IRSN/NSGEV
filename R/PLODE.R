##' Compute the lower and upper limits of a confidence interval on the
##' quantiles \eqn{q_{M^\star}(p)}{q_M(p)} of the maximum
##' \eqn{M^\star}{M} on a given period of time. This computation is
##' made for a number of probabilities \eqn{p}.
##' 
##' For a given probability \eqn{p}, the profile likelihood confidence
##' interval for the quantile
##' \eqn{\eta(\boldsymbol{\psi},\,p):=q_{M^\star}(p;\,\boldsymbol{\psi})}{eta(psi,
##' p) := qM(p, psi)} can be obtained by minimising and maximising
##' \eqn{\eta(\boldsymbol{\psi},\,p)}{eta(psi, p)}
##' w.r.t. \eqn{\boldsymbol{\psi}}{psi} under the constraint
##' \eqn{\ell(\boldsymbol{\psi}) \geq \ell_{\max} - \delta}{ell(psi)
##' >= ell_max - delta} where \eqn{\ell(\boldsymbol{\psi})}{ell(psi)}
##' is the log-likelihood function, and \eqn{\delta := q_{\chi^2(1)}(1
##' -\alpha) /2}{ qchisq(1 - alpha, df = 1) /2}.  By differentiating
##' the Lagrangian for this problem w.r.t. the probability
##' \eqn{p}, we get an ODE for the vector
##' \eqn{\boldsymbol{\psi}(p)}{psi(p)}. This ODE may be called a
##' \emph{tangential} ODE inasmuch as it describes the evolution of the
##' vector \eqn{\boldsymbol{\psi}(p)}{psi(p)} on the likelihood
##' surface defined in the parameter space by the equation
##' \eqn{\ell(\boldsymbol{\psi}) = \ell_{\max} - \delta}{ell(psi) =
##' ell_max - delta}.  The initial ODE can be found either by solving
##' a constrained optimisation problem or by solving another ODE, the
##' \emph{radial} ODE that describes the evolution of a point moving from
##' an origin at the ML estimate \eqn{\hat{\boldsymbol{\psi}}}{psiHat}
##' to a destination located on the surface likelihood.
##' 
##' @title Profile Likelihood Confidence Intervals for the Quantiles
##'     of the Maximum on a given Period
##' 
##' @param object An object with class \code{"TVGEV"}.
##'
##' @param date A vector that can be coerced to the class
##'     \code{"Date"} giving the beginnings of the years in the period.
##'
##' @param probIni The probability that will be used for the
##'     initialisation of the tangential ODE.
##' 
##' @param level The confidence level \eqn{1 - \alpha}{1 - alpha}.
##' 
##' @param out The type of output wanted. The choice \code{"array"}
##'     leads to an array with its dimensions corresponding to the
##'     probability \eqn{p}, the type of result (lower/upper bound,
##'     quantile) and the confidence level. With the choice
##'     \code{"data.frame"} we get an object inheriting from this
##'     class, actually a \code{quantMax.TVGEV} object for which some
##'     methods like \code{autoplot} are available.
##' 
##' @param trace Integer level of verbosity.
##'
##' @return An object with class \code{"quantMax.TVGEV"} inheriting
##'     from \code{"data.frame"}
##' 
##' @references
##' Yves Deville (2024)
##' "Profile Likelihood via Optimisation and Differential Equations" \emph{arXiv}
##' \doi{10.48550/arXiv.2404.02774}.
##'
##' @section Caution: This function is experimental.
##'
##' @importFrom utils packageVersion
##' 
##' @export
##' 
##' @examples
##' \dontrun{
##'    library(deSolve)
##'    df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##'    ## fit a TVGEV model. Only the location parameter is TV.
##'    object <- TVGEV(data = df, response = "TXMax", date = "Date",
##'                    design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'                    estim = "nloptr",
##'                    loc = ~ t1 + t1_1970)
##'    newDate <- as.Date(paste0(2025:2034, "-01-01"))
##'    st1 <- system.time(
##'       plODE <- quantMaxPLODE(object, date = newDate, level = c(0.70, 0.95))
##'    )
##'    st2 <- system.time(
##'        pl <- quantMax(object, date = newDate, conf = "proflik", trace = 1,
##'                   level = c(0.70, 0.95))
##'    )
##'    g <- autoplot(pl)
##'    g <- g + geom_line(aes(x = 1 - Prob, y = L, group = Level), data = plODE) +
##'       geom_line(aes(x = 1 - Prob, y = U, group = Level), data = plODE)
##'    g
##' }
quantMaxPLODE <- function(object, date, probIni = 0.70, level = 0.95,
                          out = c("data.frame", "array"), trace = 1) {
    
    if (!requireNamespace("deSolve", quietly = TRUE)) {
        stop("This function requires the use of the `deSolve' package")
    }

    if (utils::packageVersion("nieve") < "0.1.5") {
        stop("This function requires 'nieve >= 0.1.5'.",
             "Use the GitHub version if needed: https://github.com/yvesdeville/nieve")
    }
    
    out <- match.arg(out)
    p <- object$p
    newDate <- date
    
    indLevel <- order(level)
    level <- level[indLevel]
    fLevel <- formatLevel(level)
    nLevel <- length(level)
    delta <- qchisq(level, df = 1) / 2

    ## =========================================================================
    ## Function for the first ODE or "radial ODE" used with two
    ## initial conditions.  Start from a point close to the ML
    ## estimate and end on the surface likelihood with equation ell =
    ## ell_max - delta. This technique may be called "inflating the
    ## likelihood bubble".
    ## =========================================================================
    
    f1 <- function(time, state, parms) {
        psi <- state[1L:p]
        nu <- state[p + 1L]
        eta <- qMax.TVGEV(object, p = probIni, date = newDate,
                          psi = psi, deriv = TRUE, hessian = TRUE)
        nll <- object$negLogLikFun(psi, deriv = TRUE, hessian = TRUE,
                                   object = object)
        if (FALSE) {
            etap <- eta
            attr(etap ,"gradient") <- attr(etap ,"hessian") <- NULL
            print(c(eta = round(etap, digits = 3),
                    nll = round(nll$objective, digits = 3),
                    nu = round(nu, digits = 3)))
        }
        g <- -nll$gradient
        H <- drop(attr(eta, "hessian")) - nu * nll$hessian
        mat <- matrix(0, nrow = p + 1L, ncol = p + 1L)
        mat[1:p, 1:p] <- H
        mat[1:p, p + 1] <- g
        mat[p + 1, 1:p] <- g
        sec <- rep(0, p + 1)
        sec[p + 1] <- -1
        list(solve(mat, sec)) 
    }
    
    ## =========================================================================
    ## Function for the second ODE, or "tangential ODE", which is used
    ## with two initial conditions.  Start from a the solution of the
    ## first ODE and move on the surface likelihood with equation ell
    ## = ell_max - delta in the direction minimising or maximising
    ## eta(psi) := q_M(psi). See the "Computing Details" document for
    ## the crossed derivative w.r.t. p and psi.
    ## =========================================================================
    
    f2 <- function(time, state, param) {
        psi <- state[1L:p]
        nu <- state[p + 1L]
        eta <- qMax.TVGEV(object, p = time, date = newDate,
                          psi = psi, deriv = TRUE, hessian = TRUE)    
        nll <- object$negLogLikFun(psi, deriv = TRUE, hessian = TRUE,
                                   object = object)
        if (FALSE) {
            etap <- eta
            attr(etap ,"gradient") <- attr(etap ,"hessian") <- NULL
            print(c(time = time,
                    eta = round(etap, digits = 3),
                    nll = round(nll$objective, digits = 3),
                    nu = round(nu, digits = 3)))
        }
        FM <- pMax.TVGEV(object = object, q = eta, date = newDate,
                         psi = psi, deriv = TRUE)
        fM <- dMax.TVGEV(object = object, x = eta, date = newDate,
                         psi = psi, derx = TRUE, deriv = TRUE)
        derX <- -(-attr(fM, "derx") * attr(FM, "gradient") / fM +
                  attr(fM, "gradient")) / fM / fM
        g <- -nll$gradient
        H <- drop(attr(eta, "hessian")) - nu * nll$hessian
        mat <- matrix(0, nrow = p + 1L, ncol = p + 1L)
        mat[1:p, 1:p] <- H
        mat[1:p, p + 1] <- g
        mat[p + 1, 1:p] <- g
        sec <- rep(0, p + 1)
        sec[1:p] <- -derX
        list(solve(mat, sec)) 
    }
    
    nt1 <- 30
    times1 <- seq(from = 0, to = delta[nLevel] + 1e-7,
                  length.out = nt1 - nLevel)
    times1 <- c(delta, times1)
    o <- order(times1)
    times1 <- times1[o]
    ## 'times1[indLev]' is delta
    indLev <- order(o)[1:nLevel]
    
    nt2 <- 100
    times2 <- seq(from = probIni, to = 1 - 1e-4, length.out = nt2)
    fProb <- format(times2)
    Quant <- array(NA_real_,
                   dim = c(Prob = nt2, Lim = 3L, Level = nLevel),
                   dimnames = list(Prob = fProb,
                                   Type = c("Quant", "L", "U"),
                                   Level = fLevel))
    
    Quant[ , "Quant", ] <- rep(qMax.TVGEV(object = object, p = times2,
                                          date = date, psi = coef(object)),
                               nLevel)
    
    ## Initial values for the state vector psiDag := c(psi, nu) where 'nu'
    ## is the Lagrange multiplier
    PsiDagIni <- array(NA_real_,
                       dim = c(Param = p + 1, Lim = 2L, Level = nLevel + 1L),
                       dimnames = list(Param = c(parNames(object), "nu"),
                                       Type = c("L", "U"),
                                       Level = c("0", fLevel)))
    
    ## Find the initial conditions for the first differential equation
    eta <- qMax.TVGEV(object, p = probIni, date = newDate,
                      deriv = TRUE, hessian = TRUE)
    psiHat <- coef(object)
    nll <- object$negLogLikFun(psiHat, deriv = TRUE, hessian = TRUE,
                               object = object)

    ## Compute the initial values for the ODE
    delta1 <- 1e-4
    H0 <- nll$hessian
    h0 <- attr(eta, "gradient")
    h0a <- solve(H0) %*% t(h0)
    nuTilde0 <- sqrt(drop(h0 %*% h0a) / 2 / delta1)
    h0a <- drop(h0a)
    
    eps <- c("L" = -1, "U" = 1)

    for (LU in c("L", "U")) {
        
        if (trace) cat(sprintf("o Confidence bound: \"%s\"\n", LU))
 
        nuTilde_1 <- unname(eps[LU]) * nuTilde0
        psiTilde_1 <- psiHat + h0a / nuTilde_1
        PsiDagIni[ , LU, 1] <- c(psiTilde_1, nu = nuTilde_1)
            
        if (trace) {
            cat("   o Solving the first ODE for the probability `probIni`\n")
        }

        out1 <- deSolve::lsoda(y = PsiDagIni[ , LU, 1],
                               times = times1, func = f1)
        
        for (iLevel in 1:nLevel) {
            ind1i <- indLev[iLevel]
            PsiDagIni[ , LU, iLevel] <- out1[ind1i, c(parNames(object), "nu")]
            ## q[LU] <- qMax.TVGEV(object, date = newDate, p = probIni,
            ##                    psi = PsiDagIni[parNames(object), LU, iLevel])
        
            if (trace) {
                cat(sprintf("   o Solving the second ODE, level = %s\n",
                            fLevel[iLevel]))
            }
            ## solve the 2-nd ODE
            out2 <- deSolve::lsoda(y = PsiDagIni[ , LU, iLevel],
                                   times = times2, func = f2)
            for (it2 in 1:nt2) {
                Quant[it2, LU, iLevel] <-
                    qMax.TVGEV(object, date = newDate, p = times2[it2],
                               psi = out2[it2, parNames(object)])
            }
        }
       
    }
    
    if (out != "array") {
        L <- list()
        for (iLevel in seq_along(level)) {
            nm <- fLevel[iLevel]
            L[[nm]] <- data.frame(Prob = times2,
                                  ProbExc = 1.0 - times2,
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
