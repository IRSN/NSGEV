

##*****************************************************************************
##' Compute the GEV parameters for the marginal distributions of a
##' \code{TVGEV} model.
##'
##' @title Compute the Matrix of GEV Parameters from the Vector of
##' Model Parameters 'psi' and the Data
##'
##' @param model The \code{TVGEV} object.
##' 
##' @param psi A vector of parameters for \code{model}. By default the
##' vector of estimates in \code{psi} will be used.
##'
##' @param date A vector that can be coerced to the \code{"Date"}
##' class. Each element will correspond to a row of the matrix of GEV
##' parameters.
##' 
##' @param deriv Logical. If \code{TRUE} the returned object has an
##' attribute names \code{"gradient"} which is the gradient
##' array. 
##'
##' @param checkNames Logical. If \code{TRUE} it will be checked that
##' the names of \code{psi} match the \code{parNames} element of
##' \code{model}.
##'
##' @param \dots Not used yet.
##'
##' @return A matrix with \code{length(date)} rows and \code{3} colums.
##' The columns contain the location, eht scale and the shape GEV parameters
##' in that order.
##'
##' @seealso The \code{\link{GEV}} for the GEV probability functions.
##'
##' @note The \code{deriv} formal argument is mainly here to maintain
##' a compatibility with the \code{NSGEV} method. However, the
##' derivatives w.r.t. the parameter vector \code{psi} are obtained
##' simply by using the \code{X} elements of the object \code{model}.
##' 
##' 
psi2theta.TVGEV <- function(model, psi = NULL, date = NULL,
                            deriv = FALSE, checkNames = TRUE,
                            ...) {

    ## psi2theta can be used to extract the element 'theta'
    if (is.null(psi)) {
        psi <- model$estimate
    }
    
    p <- length(psi)
    parNames.GEV <- c("loc", "scale", "shape")
    
    if (checkNames) {
        if (!setequal(names(psi), parNames(model))) {
            stop("'psi' must have named elements with ",
                 "the names 'parNames(model)'")
        }
    } else {
        names(psi) <- model$parNames
    }

    ## is date is NULL, we use the design matrices in the model
    if (is.null(date)) {
        X <- model$X
        n <- model$n
        ind <- model$ind
        fDate <- model$fDate
    } else {
        n <- length(date)
        if (!all(model$isCst)) {
            L <- modelMatrices.TVGEV(model, date = date)
            X <- L$X
        } else X <- NULL
        fDate <- format(date)
    }
    
    theta <- array(NA, dim = c(n, 3L),
                       dimnames = list(rownames(X[["loc"]]), parNames.GEV))
    
    for (nm in parNames.GEV) {
        if (!model$isCst[nm]) {
            theta[ , nm] <- X[[nm]] %*% psi[model$ind[[nm]]]
        } else {
            theta[ , nm] <- psi[model$ind[[nm]]]
        }
    }
        

    if (deriv) {
        
        jac <- array(0.0, dim = c(n, model$p, 3L),
                     dimnames = list(fDate, names(psi), parNames.GEV))
        for(nm in parNames.GEV) {
            if (!model$isCst[nm]) {
                jac[ , model$ind[[nm]], nm] <- X[[nm]]
            } else {
                jac[ , model$ind[[nm]], nm] <- 1.0
            }
        }
        attr(theta, "gradient") <- jac
    }
    
    theta
    
}

##*****************************************************************************
##' Initial Parameter Estimates for a \code{TVGEV} model.
##'
##' The model parameters that are inverse-linked to the location
##' parameter are found by using a linear model with a constant error
##' variance. The initial value of the shape parameter(s) being set to
##' zero, which correspond to Gumbel margins, we know that the
##' expectation is approximately \eqn{\mu + 0.577 \sigma} while
##' the standard deviation is \eqn{\sigma \sqrt{\pi} / 6}
##' 
##' @title Initial Parameter Estimates
##' 
##' @param object A (possibly incomplete) \code{TVGEV} object. 
##'
##' @param y The response, by defaut the response attached to the
##' model.
##'
##' @return A vector of initial values.
##' 
##' 
##' 
parIni.TVGEV <- function(object, y = NULL) {
    
    if (is.null(y)) {
        y <- object$data[ , object$response]
    } else {
        if (length(y) != object$n) {
            stop("the length of 'y' must be equal to object$n")
        }
    }

    if (object$isCst["loc"]) {
        sigma <- sqrt(6) * sd(y, na.rm = TRUE) / pi
        mu <- mean(y, na.rm = TRUE) - 0.578 * sigma
        return(c("loc" = mu, "scale" = sigma, "shape" = 0.0))
    }
    
    psi0 <- list("loc" = numeric(0), "scale" = numeric(0), "shape" = numeric(0))
    psi0[["shape"]] <- rep(0.0, object$pp["shape"])
    
    dfAll <- object$dfAll
    dfAll[ , object$response] <- y
    fm <- object$loc
    fm <- update.formula(fm, as.formula(paste(object$response, " ~ .")))
    fit_loc <- lm(formula = fm, data = dfAll)
    co <- coef(fit_loc)
    
    ## the residual sd is sigma * pi / sqrt(6)
    ## XXX this works when the scale formula has an intercept
    
    psi0[["scale"]] <- rep(0.0, object$pp["scale"])
    psi0[["scale"]][1L] <- sqrt(6) * summary(fit_loc)$sigma / pi
    
    if (attr(object$terms[["loc"]], "intercept")) {
        co[1] <- co[1] - psi0[["scale"]][1] * 0.578
    } 

    psi0[["loc"]] <- co
    psi0 <- unlist(psi0)
    
    names(psi0) <- object$parNames
    psi0

}
    
##*****************************************************************************
##' Maximum Likelihood Estimation of a \code{TVGEV} model.
##'
##' @title Maximum Likelihood Estimation of a \code{TVGEV} Model 
##' 
##' @param object A (possibly incomplete) \code{TVGEV} object.
##' 
##' @param y A numeric vector giving the response to be used.
##'
##' @param psi0 Numeric vector of initial values for the parameters.
##'
##' @param estim Character giving the optimisation to be used.
##'
##' @param parTrack \code{Logical}. If \code{TRUE}, all parameters at
##' which an evaluation of the log-Likelihood will be stored and
##' returned.
##'
##' @return A list with elements that can be copied into those
##' of a \code{TVGEV} object.
##'
##' @author Yves Deville
##'
##' @section Caution: For now it is assumed that the shape parameter
##' is constant. This assumtion is used to set the bounds on the shape
##' parameter: the lower bound is \code{-0.9} and the upper bound is
##' \code{2.0}.
##' 
MLE.TVGEV <- function(object,
                      y = NULL,
                      psi0 = NULL,
                      estim = c("optim", "nloptr"),
                      parTrack = FALSE) {

    parNames.GEV <- c("loc", "scale", "shape")
     
    estim <- match.arg(estim)
    res <- list()
    
    ## =======================================================================
    ## Build response.
    ## =======================================================================

    if (is.null(y)) {
        y <- object$data[ , object$response]
    } else {
        if (length(y) != object$n) {
            stop("the length of 'y' must be equal to object$n")
        }
    }

    ## Note that we now add or modify  the 'indVal' element in 'object'
    ## and this apears in the negLogLikFun
    
    res$indVal <- object$indVal <- !is.na(y)
    yBak <- y
    y <- y[res$indVal]
    res$nobs <- sum(res$indVal)

    ## =======================================================================
    ## When parTrack is TRUE, the value of 'psi' is catched before
    ## computing the log-likelihood.
    ## =======================================================================

    if (parTrack) {
        trackEnv <- new.env()
        trackEnv$psi <- numeric(0)
    }
    
    negLogLikFun <- function(psi, deriv = FALSE, object) {
        
        ## redefine 'parNames.GEV' because this closure is returned by
        ## the estimation function.  It could be faster NOT to use
        ## names but rather integer indices.
        
        parNames.GEV <- c("loc", "scale", "shape")
        theta <- array(NA, dim = c(object$n, 3L),
                       dimnames = list(NULL, parNames.GEV))
      
        for (nm in parNames.GEV) {
            if (!object$isCst[nm]) {
                theta[ , nm] <- object$X[[nm]] %*% psi[object$ind[[nm]]]
            } else {
                theta[ , nm] <- psi[object$ind[[nm]]]
            }
        }
        
        logL <-
            dGEV(y,
                 loc = theta[object$indVal, "loc"],
                 scale = theta[object$indVal, "scale"],
                 shape = theta[object$indVal, "shape"], log = TRUE,
                 deriv = deriv)
        
        nL <- -sum(logL)
      
        if (is.na(nL)) nL <- NaN
       
        if (deriv) {
            
            grad1 <- -attr(logL, "gradient")
            gradnL <- rep(0.0, object$p)
            ## names(gradnL) <- object$parNames
            
            for (nm in parNames.GEV) {
                if (!object$isCst[nm]) {
                    gradnL[object$ind[[nm]]] <- 
                        crossprod(grad1[ , nm, drop = FALSE],
                                  object$X[[nm]][object$indVal, , drop = FALSE])
                } else {
                     gradnL[object$ind[[nm]]] <- sum(grad1[ , nm, drop = FALSE])
                }
            }
         
            return(list("objective" = nL , "gradient" = gradnL))

        }
       
        nL
        
    }
    
    res$negLogLikFun <- negLogLikFun
    
    psi0 <- res$psi0 <- parIni.TVGEV(object = object, y = yBak)
    
    if (parTrack) {
        negLogLikFun1 <- function(psi, deriv = FALSE, object)  {
            trackEnv$psi  <- c(trackEnv$psi, psi)
            negLogLikFun(psi = psi, deriv = deriv, object = object)
        }
    } else {
        negLogLikFun1 <- negLogLikFun
    }
        
    cvg <- TRUE
    if (estim == "optim") {

        res$fit <- try(optim(par = psi0,
                             fn = negLogLikFun1,
                             deriv = FALSE,
                             control = list(maxit = 1000),
                             ## hessian = TRUE,
                             object = object))

        if (!inherits(res$fit, "try-error")) {
            if (res$fit$convergence == 0) {
                estimate <- res$fit$par
                names(estimate) <- object$parNames
                res$estimate <- estimate
                res$negLogLik <- res$fit$value
            } else {
                cvg <- FALSE
            }
        }
            
    } else if (estim == "nloptr") {
        
        opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                     "xtol_rel" = 1.0e-8,
                     "xtol_abs" = 1.0e-8,
                     "ftol_abs" = 1e-5,
                     "maxeval" = 1000, "print_level" = 0,
                     "check_derivatives" = FALSE)

        ## XXX caution! this works when the shape is constant only!!!
        p <- object$p
        lb <- rep(-Inf, p)
        lb[p] <- -0.9
        ub <- rep(Inf, p)
        ub[p] <- 2.0
        
        res$fit <- try(nloptr(x0 = psi0,
                              eval_f = negLogLikFun,
                              lb = lb,
                              ub = ub,
                              deriv = TRUE,
                              opts = opts,
                              object = object))
        
        if (!inherits(res$fit, "try-error")) {
            if (res$fit$status > 0 ) {
                estimate <- res$fit$solution
                names(estimate) <- object$parNames
                res$estimate <- estimate
                res$negLogLik <- res$fit$objective
            } else {
                cvg <- FALSE
            }
        }
        
    }

    if (!cvg) {
        warning("covergence not reached in optimisation")
        estimate <- rep(NA, object$p)
        names(estimate) <- object$parNames
        res$negLogLik <- NA
        res$estimate <- estimate
        res$logLik <- NA
    } else {
        
        res$logLik <- -res$negLogLik
        psiHat <- res$estimate
        ## compute theta 
        res$theta <- psi2theta(object, psi = psiHat, checkNames = FALSE)
        
        ## compute Hessian. 
        ## res$hessian <- hessian(func = negLogLikFun,
        ##                       x = res$estimate, deriv = FALSE,
        ##                       object = object)

        res$hessian <- optimHess(par = res$estimate,
                                 fn = negLogLikFun,
                                 deriv = FALSE,
                                 object = object)
        
        vcov <- try(solve(res$hessian), silent = TRUE)
        
        if (!inherits(vcov, "try-error")) {
            rownames(vcov) <- colnames(vcov) <- object$parNames
            res$vcov <- vcov
            res$sd <- sqrt(diag(vcov))
        }
    }
    
    if (parTrack) {
        tpsi <-  matrix(trackEnv$psi, ncol = object$p,
                                byrow = TRUE)
        colnames(tpsi) <- object$parNames
        res$tracked <-
            list(psi = tpsi,
                 negLogLik = apply(tpsi, 1, negLogLikFun, deriv = FALSE,
                     object = object))
    }

    res

}

##*****************************************************************************
##' Bootstrap for a \code{TVGEV} model.
##'
##' The parametric bootstrap (\code{type = "param"}) is as follows: for
##' each bootstrap sample, a response vector
##' \eqn{\mathbf{y}^\star}{y*} is drawn from the estimated GEV
##' distribution as described in \code{object}; a MLE is performed to
##' find the bootstraped parameter vector
##' \eqn{\boldsymbol{\psi}^\star}{\psi*}. For the nonparametric
##' bootstrap (\code{type = "NP"}), the generalised residuals of the
##' object as returned by \code{\link{residuals.TVGEV}} are resampled
##' and are back-transformed to simulated observations; a MLE is then
##' performed to find the bootstraped parameter vector
##' \eqn{\boldsymbol{\psi}^\star}{\psi*} as in the parametric case.
##' 
##'
##' @title  Bootstrap for a \code{TVGEV} Model
##' 
##' @param object A \code{TVGEV} object.
##' 
##' @param R Target number of bootstrap samples.
##'
##' @param type Character. The values \code{"param"} and \code{"NP"}
##' can be used to chose between parametric and nonparametric bootstrap,
##' see \bold{Details}.
##' 
##' @param estim Argument passed to \code{MLE}.
##'
##' @param parallel Logical. If \code{TRUE} the loop over the bootstrap
##' sample is performed using \code{foreach} with \code{\%do par\%}.
##'
##' @param ... Further arguments passed to \code{MLE}.
##'
##' @return A matrix containing the bootstraped parameter vectors,
##' with one row for each bootstrap sample.
##'
##' @note The simulated response does not cover the observations where
##' the original response was \code{NA}.
##' 
##' @section Caution: For some bootstrap samples, the MLE may fail, in
##' which case the coefficients will be ignored. So the true number
##' of sampled coefficient vecors  will generally be smaller than the
##' target number as given in \code{R}.
##'
##' @examples
##'
##' example(TVGEV)
##' bs <- bs.TVGEV(res2, R = 50, estim = "nloptr")
##' 
##' \dontrun{
##'    library(parallel)
##'    library(doParallel)
##'    nc <- detectCores()
##'    cl <- makeCluster(nc)
##'    registerDoParallel(cl)
##'
##'    ## findings: with 'nloptr', less than 1% of the optimisations
##'    ## diverges, while more than 10% diverged with 'optim'.
##'    te <- system.time(bsp <- bs.TVGEV(res2, R = 5000, estim = "nloptr",
##'                                      parallel = TRUE))
##'    stopCluster(cl)
##' }
##'
##' 
bs.TVGEV <- function(object,
                     R = 100,
                     type = c("param", "NP"),
                     estim = "optim",
                     parallel = FALSE,
                     ...) {

    type <- match.arg(type)
    Psi <- array(NA, dim = c(R, object$p),
                 dimnames = list(NULL, object$parNames))
    
    nlL <- rep(NA, R)

    ## ==================================================================
    ## CAUTION: 'y' should be NA where the original response was NA
    ## otherwise the uncertainty on the estimated parameters will be
    ## too small.
    ## =================================================================

    iv <- object$indVal
    nv <- sum(iv)
    y <- array(NA, dim = c(object$n, R),
               dimnames = list(NULL, paste("sim", 1:R, sep = "_")))
    
    if (type == "param") {
        
        y <- array(NA, dim = c(object$n, R),
               dimnames = list(NULL, paste("sim", 1:R, sep = "_")))
        
        y[iv, ]  <- rGEV(n = R,
                         loc = object$theta[iv, "loc"],
                         scale = object$theta[iv, "scale"],
                         shape = object$theta[iv, "shape"])
        
    } else if (type == "NP") {
       
        e <- array(unname(resid(object, type = "unif")[iv]),
                   dim = c(nv, R))

        samp <- function(x) {
            qGEV(p = sample(x, size = nv, replace = TRUE),
                 loc = object$theta[iv, "loc"],
                 scale = object$theta[iv, "scale"],
                 shape = object$theta[iv, "shape"])
        }
        
        y[iv, ] <- apply(X = e, MARGIN = 2, FUN = samp) 
    }
    
    if (parallel) {
        
        if (requireNamespace("foreach", quietly = TRUE)) {
            Psi <-
                foreach::"%dopar%"(foreach::foreach(b = 1:ncol(y),
                                                    .export = c("MLE.TVGEV", "parIni.TVGEV"),
                                                    .packages = c("nloptr", "NSGEV", "numDeriv"),
                                                    .combine = "rbind"), {
                                       res <- try(MLE.TVGEV(object, y = y[ , b], estim = estim, ...),
                                                  silent = TRUE)
                                       if (!inherits(res, "try-error")) {
                                           est <- res$estimate
                                       } else est <- NULL
                                       est
                                   })
            
            nlL <- rep(NA, ncol(y))
        } else {
            stop("the package 'foreach' could not be used")
        }
            
    } else {
        
        for (b in 1:R) {
            res <- try(MLE.TVGEV(object, y = y[ , b], estim = estim, ...),
                       silent = TRUE)
            if (!inherits(res, "try-error")) {
                Psi[b, ] <- res$estimate
                nlL[b] <- res$negLogLik
            }
        }

        ## XXX return the negative log-likelihood as well?
    }
    
    ind <- apply(Psi, 1, function(x) !any(is.na(x))) 
    Psi <- Psi[ind, , drop = FALSE]
    nlL <- nlL[ind]
    
    list(estimate = Psi,
         negLogLik = nlL,
         R = R,
         optim = optim,
         type = type)
    
}

## *************************************************************************
##' Build model matrices.
##'
##' These matrices are needed when the model is build (and usually is
##' estimated) or when a prediction is required.
##' 
##' @title Build Model Matrices
##'
##' @param object An object with S3 class \code{"TVGEV"}.
##'
##' @param date An object with class \code{"Date"} or an object that
##' can be coerced to this class.
##'
##' @return A list with the follwing elements
##'
##' \item{dfAll}{
##'
##' A data frame with all the variables required.
##'
##' }
##' \item{X}{
##'
##' A list with three elements named \code{"loc"}, \code{"scale"} and
##' \code{"shape"}. Each element is the corresponding model matrix in
##' the \code{lm} meaning. Each matrix can be multiplied by the
##' corresponding coefficients to get the prediction for the
##' corresponding parameter.
##'
##' }
##'
##' @seealso The function \code{\link[stats]{model.matrix}} used by
##' \code{\link[stats]{lm}}.
##' 
modelMatrices.TVGEV  <- function(object, date = NULL) {
    
    
    if (is.null(date)) {
        return(list(dfAll = object$dfAll, X = object$X))
    }

    date <- as.Date(date)
    data <- data.frame(date)
    colnames(data) <- object$date
    parNames.GEV <- c("loc", "scale", "shape")

    ## Caution: here, 'dfAll' does not embed the response
    if (!is.null(object$design)) {
        if (is.call(object$design)) {
            dfAll <- as.data.frame(eval(object[["design"]], envir = data))
        } else {
            dfAll <- object$design
        }
    } else dfAll <- object$data

    fDate <- format(date)
    rownames(dfAll) <- fDate
    
    X <- list()
    
    for (nm in parNames.GEV) {
        X[[nm]] <- model.matrix(object$terms[[nm]], data = dfAll)
        rownames(X[[nm]]) <- fDate
    }

    list(dfAll = dfAll, X = X)

}

## *************************************************************************
##' Time-varying GEV model
##'
##' This kind of model describe \emph{independent} observations having
##' a Generalized Extreme Value (GEV) distribution depending on time.
##' The three GEV parameters (\code{"location"}, \code{"scale"} and
##' \code{"shape"}) can depend on a date variable used as a covariate
##' in a linear-model style of dependence. The covariates of the model
##' are functions of the date as provided in \code{breaksX},
##' \code{polynomX}, \code{trigonX}. No other covariates can be used
##' unless the use of \code{predict} or \code{RL} will not be possible.
##' 
##' @title Time-Varying GEV Model.
##' 
##' @param data A data frame containing at least the two required
##' variables with their names given in \code{date} and
##' \code{response}.
##'
##' @param date Character. Name of the date variable in \code{data}.
##'
##' @param response Character. Name of the response variable in
##' \code{data}.
##'
##' @param design A call to a function creating a data frame with all
##' the variables required in the formulas \code{loc}, \code{scale}
##' and \code{shape}. While the content of \code{data} is mainly used
##' for the creation of the \code{TVGEV} object, \code{design} will be
##' used in prediction or computation of Return Levels.
##' 
##' @param loc A formula linking linearly the GEV location parameter
##' to the columns in the data frame created with \code{design}. The
##' Rogers-Wilkinson style of formula is used here (as in \code{lm}):
##' only the variables are given, the parameter names are created in
##' relation to the variables names.
##'
##' @param scale A formula for the GEV scale parameter.
##' 
##' @param shape A formula for the GEV shape parameter.
##'
##' @param psi0 A numeric vector of initial values for the parameters.
##' By default a vector is computed using a linear model.
##' 
##' @param estim Character giving the package or function to be used
##' for the likelihood maximisation.
##'
##' @param parTrack Logical. If \code{TRUE} the value of the parameter
##' vector is tracked during the optimisation (i.e. at each call of
##' the objective) and all the values are returned as a matrix
##' \code{"psi"} in a list named \code{"tracked"} which contains as
##' well the value of the objective.
##' 
##' @param trace Integer level of verbosity.
##'
##' @section Caution: The call passed in the \code{design} formal will
##' be reused in prediction. Passing a merely data frame would make
##' prediction impossible because there would then be no way to guess
##' how the covariates of the model were created from the date.
##' 
##' @return An object with class \code{"TVGEV"}, a list. Among the
##' list elements
##'
##' \item{call}{The call. This will be required.}
##' 
##' \item{data, design, loc, scale, shape}{Copies of the inputs.}
##'
##' \item{response, date}{The \emph{name} of the response and of the
##' date columns in \code{data}.}
##'
##' \item{terms}{A list with the terms for the three GEV parameters.}
##' 
##' \item{opt}{The results of the optimisation, if relevant.}
##'
##' \item{estimate}{The numeric vector of estimates for the vector
##' \code{psi}.}
##' 
##' \item{theta}{A matrix with three columns containing the GEV
##' parameters \code{loc}, \code{scale} and \code{shape}.}
##' 
##' @note When the response contains NA, the corresponding elements
##' are discarded in the computation of the log-likelihood. However the
##' corresponding rows exist in the matrix \code{theta}.
##' 
##' @author Yves Deville
##'
##' @examples
##'
##' ## transform a numeric year into a date
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' df0 <- subset(df, !is.na(TXMax))
##' 
##' ## fit a TVGEV model. Only the location parameter is TV.
##' t1 <- system.time(
##'     res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'                   design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'                   loc = ~ t1 + t1_1970))
##' 
##' ## The same using "nloptr" optimisation.
##' t2 <- system.time(
##'     res2 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'                   design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'                   loc = ~ t1 + t1_1970,
##'                   estim = "nloptr",
##'                   parTrack = TRUE))
##' 
##' ## use extRemes::fevd the required variables need to be added to the data frame
##' ## passed as 'data' argument
##' t0 <- system.time({
##'    df0.evd <- cbind(df0, breaksX(date = df0$Date, breaks = "1970-01-01",
##'                     degree = 1));
##'    res0 <- fevd(x = df0.evd$TXMax, data = df0.evd, loc = ~ t1 + t1_1970)
##'  })
##'
##' ## compare estimate and negative log-liks
##' cbind("fevd" = res0$results$par,
##'       "TVGEV_optim" = res1$estimate,
##'       "TVGEV_nloptr" = res2$estimate)
##' cbind("fevd" = res0$results$value,
##'       "VGEV_optim" = res1$negLogLik,
##'       "TVGEV_nloptr" = res2$negLogLik)
##' 
##' ## ====================================================================
##' ## use a loop on plausible break years. The fitted models
##' ## are stored within a list
##' ## ====================================================================
##'
##' \dontrun{
##' 
##'     yearBreaks <- c(1940, 1950, 1955, 1960:2000, 2005, 2010)
##'     res <- list()
##' 
##'     for (ib in seq_along(yearBreaks)) {
##'         d <- sprintf("%4d-01-01", yearBreaks[[ib]])
##'         floc <- as.formula(sprintf("~ t1 + t1_%4d", yearBreaks[[ib]]))
##'         res[[d]] <- TVGEV(data = df, response = "TXMax", date = "Date",
##'         design = breaksX(date = Date, breaks = d, degree = 1),
##'         loc = floc)
##'     }
##'
##'     ## [continuing...] ]find the model with maximum likelihood, and plot
##'     ## something like a profile likelihood for the break date considered
##'     ## as a new parameter. However, the model is not differentiable w.r.t.
##'     ## the break! 
##' 
##'     ll <- sapply(res, logLik)
##'     plot(yearBreaks, ll, type = "o", pch = 21, col = "orangered",
##'          lwd = 2, bg = "gold", xlab = "break", ylab = "log-lik")
##'     grid()
##'     iMax <- which.max(ll)
##'     abline(v = yearBreaks[iMax])
##'     abline(h = ll[iMax] - c(0, qchisq(0.95, df = 1) /2),
##'            col = "SpringGreen3", lwd = 2)
##'
##' }

TVGEV <- function(data,
                  date,
                  response,
                  design = NULL,
                  loc = ~ 1,
                  scale = ~ 1,
                  shape = ~ 1,
                  psi0 = NULL,
                  estim = c("optim", "nloptr", "none"),
                  parTrack = FALSE,
                  trace = 0) {

    estim <- match.arg(estim)
    parNames.GEV <- c("loc", "scale", "shape")
    symNames.GEV <- c("loc" = "mu", "scale" = "sigma", "shape" = "xi")
    
    data <- as.data.frame(data)
    n <- nrow(data)
    
    if (!(response %in% names(data))) {
        stop("'response' is not a colname of 'data'")
    }

    ## =======================================================================
    ## Start building the model
    ## =======================================================================
    
    mc <- match.call()
    tv <- list(call = mc, data = data)

    ## =======================================================================
    ## Build a data frame 'dfAll' containing ALL the required variables.
    ## =======================================================================

    tv[["design"]] <- mc[["design"]]

    if (!is.null(tv$design)) {
        if (is.call(tv$design)) {
            dfAll <- as.data.frame(eval(tv$design, envir = data))
            dfAll <- data.frame(data[ , response], dfAll)
            colnames(dfAll) <- c(response, colnames(dfAll)[-1])
        } else {
            dfAll <- tv$design
        }
    } else dfAll <- tv$data

    tv$fDate <- format(data[ , date])
    rownames(dfAll) <- tv$fDate
    
    ## =======================================================================
    ## for each GEV parameter, get the formula given under that name
    ## and build a data frame with the required variables stored in
    ## the list 'X'. Note that for a (default) formula '~ 1' a data
    ## frame with zero columns is returned by 'model.frame' but a
    ## matrix with ones is returned by 'model.matrix'.
    ## =======================================================================

    isCst <- c("loc" = FALSE, "scale" = FALSE, shape = FALSE)
    pp <- c("loc" = 0L, scale = 0L, shape = 0L)
    ppPrev <- 0L
    X <- parNames <- psi0 <- trm <- ind <- list()
    
    for (nm in parNames.GEV) {

        fm <- tv[[nm]] <- get(nm)
        
        trm[[nm]] <- terms(fm)
        pn <- attr(trm[[nm]], "term.labels")
        if (length(pn) == 0) isCst[nm] <- TRUE
    
        if (attr(trm[[nm]], "intercept")) pn <- c("0", pn)
        pn <- paste(symNames.GEV[nm], pn, sep = "_")
        parNames[[nm]] <- pn
        pp[nm] <- length(pn)
        
        psi0[[nm]] <- rep(NA, length(pn))
        
        X[[nm]] <- model.matrix(trm[[nm]], data = dfAll)
        ind[[nm]] <- ppPrev + (1L:pp[[nm]])
        ppPrev <- ppPrev + pp[[nm]]
    }
    
    parNames <- unlist(parNames)
    p <- length(parNames)
    
    tv <- c(tv,
            list(response = response, date = date,
                 dfAll = dfAll, terms = trm, X = X,
                 isCst = isCst,
                 n = n, p = p, pp = pp, ind = ind, parNames = parNames,
                 df = p))
  
    ## tv$theta <- attr(negLogLik(, "theta")
    class(tv) <- "TVGEV"

    res <- MLE.TVGEV(object = tv,
                     y = NULL,
                     psi0 = psi0,
                     estim = estim,
                     parTrack = parTrack) 

    ## copy
    for (nm1 in names(res)) {
        tv[[nm1]] <- res[[nm1]]
    }
    
    tv$theta <- psi2theta(tv, psi = tv$estimate)
    rownames(tv$theta) <- tv$fDate

    tv
    
}

## *************************************************************************
##' Random simulation from a \code{TVGEV} object.
##' 
##' @title Simulate from a \code{TVGEV} Object
##'
##' @param object An object of class \code{"TVGEV"} representing a
##' Time-Varying GEV model.
##' 
##' @param nsim Number of simulated 
##' 
##' @param seed Not used yet.
##'
##' @param newdate A vector with class \code{"Date"} or that can be
##' coerced to this class. The default \code{NULL} leads to using the
##' date used in \code{object}.
##'
##' @param psi A vector
##'
##' @param ... Not used yet.
##'
##' @return A matrix with \code{nsim} columns and one row by date.
##' This matrix is given a special S3 class \code{"simulate.TVGEV"},
##' mainly in order to facilitate plotting.
##'
##' @examples
##'
##' example(TVGEV)
##' sim <- simulate(res2, nsim = 200)
##' plot(sim)
##' 
simulate.TVGEV <- function (object, nsim = 1, seed = NULL,
                            newdate = NULL, psi = NULL,  ...) {
    
    theta <- psi2theta(model = object, psi = psi, date = newdate)
    
    sim <- rGEV(nsim, loc = theta[, 1L], scale = theta[, 2L], 
                shape = theta[, 3L])
    
    dimnames(sim) <- list(rownames(theta),
                          paste("sim", 1L:nsim,  sep = ""))

    if (is.null(newdate)) newdate <- object$data[ , object$date]
    
    attr(sim, "date") <- newdate
    
    class(sim) <- "simulate.TVGEV"
    
    sim

}

## *************************************************************************
##' Plot Paths Simultated from a \code{TVGEV} object
##'
##'
##' @title Plot Paths Simultated from a \code{TVGEV} object
##' 
##' @param x A \code{TVGEV} object
##' 
##' @param y Not used.
##' 
##' @param col Color to be used.
##'
##' @param alpha Opacity level.
##'
##' @param ... Other arguments to be passed to \code{plot}.
##'
##' @return Nothing.
##'
##' @seealso \code{\link{simulate.TVGEV}}.
##'
##' @section Caution: This function will soon be removed. The
##' \code{simulate} method will return a block time series
##' with class \code{"bts"} instead of an object with a specific
##' class \code{"simulate.TVGEV"}. This is intended to avoid
##' an unnecessary proliferation of classes.
##' 
plot.simulate.TVGEV <- function(x, y, col = "gray",
                                alpha = NULL, ...) {

    if (!missing(y) && !is.null(y)) {
        warning("'y' formal not used by this method")
    }
    
    date <- attr(x, "date")

    plot(date, x[ , 1], type = "n",
         ylim = range(x), 
         xlab = "date",
         ylab = "simulated data",
         ...)

    if (is.null(alpha)) {
        alpha <- approx(x = c(1, 10, 100, 1000, 10000),
                        y = c(1.0, 0.8, 0.5, 0.2, 0.1),
                        xout = ncol(x))$y
    }
    
    matlines(date, x, type = "o",
             pch = 16, cex = 0.6,
             col = translude(col, alpha),
             ...)

}

## *************************************************************************
## Extract the vector of coefficients from a \code{TVGEV} object.
##
## @title Coefficients of a \code{TVGEV} object
##
## @param object A \code{TVGEV} object.
## 
## @param type Character. When \code{"psi"}, the vector of model
## parameters is returned. When instead \code{type} is \code{"theta"},
## the matrix of GEV parameters is returned, with one row by block
## (or observation) and one column for each of the GEV parameters
## \code{"loc"}, \code{"scale"} and \code{"shape"}.
##
## @param ... Not used yet.
## 
## @return Vector \eqn{\mathbf{\psi}}{\psi} of coefficients, or
## matrix with the GEV parameters \eqn{\mathbf{\theta}_i}{\theta_i}
## as its rows.
## 
coef.TVGEV <- function(object, type = c("psi", "theta"),  ...) {
    type <- match.arg(type)
    if (type == "psi") return(object$estimate)
    else {
        co <- object$theta
        rownames(co) <- attr(co, "date") <- object$fDate
        attr(co, "collabels") <- c("loc", "scale", "shape")
        attr(co, "label") <- "GEV parameters"
        class(co) <- c("bts", "matrix")
        return(co)
    }
}

vcov.TVGEV <- function(object, ...) {
    object$vcov
}

logLik.TVGEV <- function(object, ...) {
    res <- object$logLik
    attr(res, "df") <- object$df
    attr(res, "nobs") <- object$nobs
    res
}

summary.TVGEV <- function(object, ...) {
    res <- object
    class(res) <- "summary.TVGEV"
    res
}

print.summary.TVGEV <- function(x, ...) {
    
    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Coefficients:\n")
    grp <- rep(c("loc", "scale", "shape"), times = x$pp)
    print(cbind("Estimate" = x$estimate,
                "Std. Error" = x$sd))
    cat("\n")
    
    cat(sprintf("Negative log-likelihood:\n%7.3f\n\n", x$negLogLik))

}
   
print.TVGEV <- function(x, ...) {
    print(summary(x))
}
