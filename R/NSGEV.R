##*****************************************************************************
##' Compute the matrix of GEV parameters from the vector of model
##' parameters 'psi' and the Data.
##'
##' @title Compute the Matrix of GEV Parameters from the Vector of
##' Model Parameters 'psi' and the Data
##' 
##' @param psi The vector of model parameters.
##'
##' @param model The GEV model: an object with class \code{"NSGEV"}.
##' 
##' @param data The data frame of covariates. By default, the data
##' frame attached to the model is used.
##'
##' @param deriv Logical. If true, the jacobian is computed and
##' returned as an attribute of the results. It is an array with
##' dimension \code{c(n, p, 3L)}
##'
##' @param checkNames Logical. If \code{TRUE}, the names of the vector
##' \code{psi} are checked against the parnames of the model.
##'
##' @param ... Not used yet.
##' 
##' @return A vector of parameters for the \code{NSGEV} model.
##' 
##' @author Yves Deville
##'
##' @method psi2theta NSGEV
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##' fit <- NSGEV(formulas = list("loc" = ~ alpha + beta * t, "scale" = ~ delta, "shape" = ~ xi),
##' data = df)
##' df.new <- data.frame(t = 11:20)
##' psi <- c("alpha" = 1, "beta" = 0.01, "delta" = 0.6, "xi" = 0.06)
##' theta.new <- psi2theta(model = fit, psi = psi, data = df.new)
##' matplot(df.new$t, theta.new, type = "b")
psi2theta.NSGEV <- function(model, psi, data = NULL, deriv = TRUE,
                            checkNames = TRUE, ...) {
   
   ## add one NAMES column for each parameter psi_1
   p <- length(psi)
   if (is.null(data)) {
      data <- model$data
   } else {
      ## checkNSGEVData(data, model)
   }
   n <- nrow(data)
   
   if (checkNames) {
      if (!setequal(names(psi), parNames(model))) {
         stop("'psi' must have named elements with ",
              "the names 'parNames(model)'")
      }
   } else {
      names(psi) <- model$parNames
   }
   
   parNames.gev <- c("loc", "scale", "shape")
   
   if (any(parNames.gev %in% names(data))) {
      stop("'data' can not have a column with its name ",
           "in ", parNames.gev)
   }
   
   env <- new.env(hash = TRUE)
   
   for(nm in names(data)) {
      assign(nm, data[[nm]], envir = env)
   }
   for(nm in names(psi)) {
      assign(nm, rep(psi[nm], n), envir = env)
   }
   for(nm in parNames.gev) {
      assign(nm, value = eval(model$der[[nm]], envir = env), envir = env)
   }
   
   res <- cbind("loc" = get("loc", envir = env),
                "scale" = get("scale", envir = env),
                "shape" = get("shape", envir = env))
   
   rownames(res) <- rownames(data)
   
   if (deriv) {
      jac <- array(NA, dim = c(n, p, 3),
                   dimnames = list(rownames(data), names(psi), parNames.gev))
      for(nm in parNames.gev) {
         mat <- get(nm, envir = env)
         jac[ , , nm] <- attr(mat, "gradient")  
      }
      attr(res, "gradient") <- jac
   }
   
   res
   
}

##*****************************************************************************
##' Find the complete vector of NSGEV parameters from a partial vector
##' and a Return Level.
##'
##' In a Non-Sationary framework, a Return Level (RL) \code{rho}
##' relates to a given value of the covariates or to a distribution of
##' the covariates. This distribution is assumed here to be given by
##' the observations in \code{data}, assuming that the corresponding
##' Return Period is equal to the number of observations understood as
##' a multiple of the block duration. When \code{type = "expect"} the
##' RL is the value \eqn{\rho} for which the random number of
##' exceedances over \eqn{\rho} has unit expectation.
##'
##' This function is a technical function to profile the likelihood.
##' One of the NSGEV model parameter, say \eqn{\psi_1}, is adjusted to
##' reach the given value RL, the other elements of \eqn{\psi} being
##' given and fixed. These later elements form a vector
##' \eqn{\boldsymbol{\psi}_{-1}}{\psi_{-1}} with length \eqn{p-1},
##' where \eqn{p} is the number of model parameters. Now for a given
##' value of \eqn{rho}, the value of the profile log-likelihood
##' \eqn{\ell(\rho)}{l(\rho)} is obtained by maximising the
##' log-likelihood w.r.t \eqn{\boldsymbol{\psi}_{-1}}{\psi_{-1}}.
##' 
##' @title Find the complete vector of NSGEV parameters from a partial
##' vector and a Return Level
##'
##' @param rho A fixed Return Level (RL).
##' 
##' @param nm1 Name of the parameter which is adjusted.
##' 
##' @param psi_m1 Vector of NSGEV parameters with its element
##' \code{nm1} removed.
##' 
##' @param model The \code{NSGEV} model.
##'
##' @param data Data frame of covariates.
##'
##' @param type The type of Return Level. 
##'
##' @return A vector \eqn{\boldsymbol{\psi}}{\psi} of NSGEV parameters.
##' 
##' @author Yves Deville
##'
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##' fit <- NSGEV(formulas = list("loc" = ~ alpha + beta * t,  "scale" = ~ delta, "shape" = ~ xi),
##'                              data = df)
##' df.new <- data.frame(t = 11:20)
##' psi <- c("alpha" = 1, "beta" = 0.01, "delta" = 0.6, "xi" = 0.06)
##' rho2psi(rho = 30, nm1 = "alpha", psi_m1 = psi[-1], model = fit, data = df.new)
##' rho2psi(rho = 40, nm1 = "alpha", psi_m1 = psi[-1], model = fit, data = df.new)
rho2psi <- function(rho, nm1, psi_m1, model, data = NULL,
                     type = "expect") {
   
   DEBUG <- FALSE
   
   if (is.null(data)) data <- model$data
   else {
       ## checkNSGEVData(data, model)
   }
   
   pnms <- parNames(model)
   p <- length(pnms)
   n <- nrow(data)
   
   ## check that 'nm1' is a coorect name
   if (!(nm1 %in% pnms)) {
      stop("'nm1' does not mach one parameter name of the model ",
           "given in 'model'")
   }
   
   ind <- match(nm1, parNames(model))
   
   if (!setequal(names(psi_m1), pnms[-ind])) {
      stop("'psi_m1' is expected to have the folowing names: ",
           pnms[-ind])
   }

   ## XXX TODO
   ## check that none of the gev parameters has a negative derivative
   ## w.r.t.  parameter 'nm1' (on the provided data).

   ## This function takes the value zero for the wanted psi_1
   gPsi <- function(psi_1) {
      psi <- rep(NA, p)
      psi[ind] <- psi_1
      psi[-ind] <- psi_m1
      names(psi) <- pnms
      theta <- psi2theta(model = model, psi = psi, data = data)
      F <- nieve::pGEV(rho, loc = theta[ , 1L], scale = theta[ , 2L],
                       shape = theta[ , 3L])
      s <- sum(F) - (n - 1.0)
      attr(s, "psi") <- psi
      if (DEBUG) {
         cat(sprintf("s = %6.2f\n\n", s))
      }
      s
   }
   
   res <- uniroot(f = gPsi, interval = c(0.01, 200))
   
   ## print(res$f.root)
   if (abs(res$f.root) > 1e-5) {
      stop("no solution found")
   }
   ## x <- res$root
   ## names(x) <- nm1
   attr(res$f.root, "psi")
}

##*****************************************************************************
##' Evaluate the negative log-likelihood of a NSGEV model.
##'
##' @title Evaluate the negative log-likelihood of a NSGEV model
##'
##' @param psi Vector of NSGEV parameters.
##'
##' @param model The NSGEV model.
##'
##' @param data A data frame with the covariates to be used. By
##' default, the data attached to \code{model} is udes.
##'
##' @param y The response as in \code{NSGEV}. By default the response
##' attached to the model is used.
##'
##' @param deriv Logical. When \code{TRUE}, the gradient
##' \eqn{-\partial \ell/\partial \boldsymbol{\psi}}{-d log L/ d \psi}
##' of the negative log-likelihood w.r.t. the parameters of the model
##' will be returned. In this case, a list will be returned.
##'
##' @param checkNames Logical. If \code{TRUE}, the names of the vector
##' \code{psi} are checked against the parnames of the model.
##' 
##' @return The negative log-likelihood (to be minimized) or a list
##' with two elements corresponding to the value of the log-likelihood
##' and its gradient.
##' 
##' @author Yves Deville
##'
##' @note The returned value of this function when \code{deriv} is \code{TRUE}
##' is designed for use with the \code{nloptr} package.
##'
##' @section Caution: The class of the result may be changed in future
##' versions.
##' 
negLogLik <- function(psi, model, data = NULL, y = NULL,
                      deriv = TRUE,
                      checkNames = TRUE) {

    if (is.null(data)) data <- model$data
    if (is.null(y)) y <- model$response
    
    theta <- psi2theta(model = model, psi = psi, data = data, deriv = deriv,
                       checkNames = checkNames) 
    n <- nrow(data)
    p <- length(psi)
    
    logL <- nieve::dGEV(y, loc = theta[, 1L], scale = theta[, 2L],
                        shape = theta[, 3L], log = TRUE, deriv = deriv)
    nl <- -sum(logL)
    if (deriv) {
        ## XXX a optimiser plus tard
        grad1 <- -attr(logL, "gradient")
        gradnl <- rep(0, p)
        names(gradnl) <- parNames(model)
        G <- attr(theta, "gradient")
        for (b in 1:n) {
            gradnl <- gradnl + grad1[b, , drop = FALSE] %*% t(G[b, , ])
        }
        return(list("objective" = nl , "gradient" = gradnl))
    } else {
        return(nl)
    }
   
}

##*****************************************************************************
##' Non-Stationary GEV.
##'
##' The model involves a vector \eqn{\boldsymbol{\psi}}{\psi} of
##' \eqn{p} parameters which are called the \emph{coefficients} of the
##' model. The GEV parameters \eqn{\mu} (loc), \eqn{\sigma} (scale)
##' and \eqn{\xi} (shape) are expressed as formulas using the
##' parameters in \eqn{\boldsymbol{\psi}}{\psi} and the model
##' variables. The three GEV parameters form a vector
##' \eqn{\boldsymbol{\theta}(\mathbf{x}) =
##' [\mu(\mathbf{x}),\,\sigma(\mathbf{x}),\,xi(\mathbf{x})]^\top}{\theta(x) = [\mu(x), \sigma(x), \xi(x)]}
##' which depends on the covariates hence varies across the
##' observations. The observations are assumed to be independent
##' conditional on the covariates.
##'
##' Note that the GEV parameters are always assumed to be given in the
##' order \code{loc}, \code{scale}, \code{shape}. The data frame \code{data}
##' can not for now use these names for its columns.
##' 
##' 
##' @title Non-Stationary GEV Model
##' 
##' @param formulas A named list with three formulas. 
##'
##' @param data The data frame to be used.
##'
##' @param response Name of the column in \code{data} that gives the
##' response. If \code{NULL}, the model will not be estimated but still
##' can have a limited use.
##'
##' @param psi Named vector of coefficients. If the model is
##' (re)-estimated, this vector is used a initial value for the
##' coefficients.
##'
##' @param est Character or logical. The default is to use the parameter
##' values given in \code{psi}.
##' 
##' @param trace Integer level of verbosity.
##'
##' @return An object with class \code{"NSGEV"}
##'
##' @author Yves Deville
##'
##' @seealso \code{\link{simulate.NSGEV}}
##'
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##'
##' ## built a model with given coefficients
##' psi <- c("alpha" = 1, "beta" = 0.01, "delta" = 0.6, "xi" = 0.06)
##' ns0 <- NSGEV(formulas = list("loc" = ~ alpha + beta * t, "scale" = ~ delta, "shape" = ~ xi),
##'              data = df, psi = psi)
##'
##' ## simulate a path
##' set.seed(1234)
##' ysim <- simulate(ns0, nsim = 1, psi = psi)
##' df2 <- cbind(df, y = ysim[ , 1L])
##' ns1 <- NSGEV(formulas = list("loc" = ~ alpha + beta * t, "scale" = ~ delta, "shape" = ~ xi),
##'              data = df2, response = "y", psi = psi, est = "optim")
##'
##' ## try an exponential link
##' ns2 <- NSGEV(formulas = list("loc" = ~ exp(alpha + beta * t), "scale" = ~ delta, "shape" = ~ xi),
##'              data = df2, response = "y", psi = psi, est = "optim")
##' 
##' ## compare the estimation with that of ismev::gev.fit
##' require(ismev)
##' ns1.ismev <- gev.fit(xdat = df2$y, ydat = as.matrix(df), mul = 1, show = FALSE)
##' rbind("NSGEV" = c(ns1$estimate, "negLogL" = ns1$negLogL),
##'        "ismev" = c(ns1.ismev$mle, "negLogL" = ns1.ismev$nllh))
##'
##' ## Try an expoential link
##' ns2.ismev <- gev.fit(xdat = df2$y, ydat = as.matrix(df), mul = 1, mulink = exp, show = FALSE)
##' rbind("NSGEV" = c(ns2$estimate, "negLogL" = ns2$negLogL),
##'        "ismev" = c(ns2.ismev$mle, "negLogL" = ns2.ismev$nllh))
NSGEV <- function(formulas,
                  data,
                  response = NULL,
                  psi = NULL,
                  ## parLower = , parUpper =
                  est = c("none", "optim", "nloptr"),
                  trace = 0) {
   
   ## XXX check the absence of circular reference
   ## among the parameters.
   
   if (!is.null(response)) {
      if (!is.character(response) ||
             !(response %in% names(data)) ||
             !is.numeric(data[[response]])) {
         stop("'response' must be the name of a numeric column ",
              "in 'data'")
      }
      y <- data[[response]]
   }
   
   data <- as.data.frame(data)
   
   parNames.gev <- c("loc", "scale", "shape")
   
   ## XXXX check names 'loc', 'scale', 'shape'
   ## in that order.
   if (!is.list(formulas) || 
          !all(sapply(formulas, is, class2 = "formula"))) {
      stop("'formulas' must be a named list ",
           "of three formulas")
   }
   
   if (any(parNames.gev %in% names(data))) {
      stop("'data' can not have a column with its name ",
           "in ", parNames.gev)
   }
   
   parNames <- list()
   allVars <- character(0)
   
   for (nm in parNames.gev) {
       av <- all.vars(formulas[[nm]])
       allVars <- union(allVars, intersect(av, names(data)))
       parNames[[nm]] <- setdiff(av, names(data))
   }
 
   parNames <- unique(unlist(parNames))
   parNames <- setdiff(parNames, parNames.gev)
   
   der <- list()
   for (nm in parNames.gev) {
      der[[nm]] <- deriv(formulas[[nm]], parNames)
   }
   
   ns <- list(formulas = formulas,
              data = data,
              der = der,
              predNames = allVars,
              p = length(parNames),
              parNames = parNames)

   class(ns) <- "NSGEV"
   
   ##==========================================================================
   ## Perform Maximum Likelihood estimation
   ##==========================================================================
   est <- match.arg(est)

   if ( (is.character(est) && est != "none") ||
       (is.logical(est) && est)) {
       estFlag <- TRUE
       if (is.null(response)) {
           stop("when 'response' is not given no estimation is done")
       }
   } else estFlag <- FALSE

   if (estFlag) {
       if (trace) cat("Maximum Likelihood Estimation\n")
       par0 <- rep(0.0, length(parNames))
       names(par0) <- parNames
       par0[names(psi)] <- psi
       
       if (trace) {
           cat("   Initial values of parameters\n")
           print(par0)
       }
       
       if (est == "nloptr") {
           
           opts <- list("algorithm" = "NLOPT_LD_LBFGS",
                        "xtol_rel" = 1.0e-8,
                        "maxeval" = 1000,
                        "print_level" = 2)
           
           fit <- nloptr(x0 = par0,
                         eval_f = negLogLik,
                         ## lb = lb
                         ## ub = c(Inf, Inf, 1.4),
                         model = ns,
                         data = data,
                         deriv = TRUE,
                         checkNames = FALSE,
                         y = y,
                         opts = opts)
           
           estimate <- fit$solution
           names(estimate) <- parNames
           ns$fit <- fit
           ns$estimate <- estimate
           ns$negLogL <- fit$objective
           
       } else if (est == "optim") {
           
           fit <- optim(par = par0, fn = negLogLik,
                        model = ns, data = data, deriv = FALSE, y = y)
           estimate <- fit$par
           names(estimate) <- parNames
           ns$fit <- fit
           ns$estimate <- estimate
           ns$negLogL <- fit$value
       }  
   } else {
       if (is.null(psi)) {
           psi <-rep(NA, ns$p)
           names(psi) <- parNames
           ns$estimate <- psi 
       } else if (!setequal(names(psi), parNames)) {
           stop("'psi' must be a named vector with suitable names")
       } else {
           ns$estimate <- psi[parNames]
           names(ns$estimate) <- parNames
       }
   }

   
   theta <- psi2theta(model = ns, psi = psi, data = ns$data, deriv = TRUE)
   thetaGrad <- attr(theta, "gradient")
   attr(theta, "gradient") <- NULL
   ns$theta <- theta
   ns$thetaGrad <- thetaGrad
   
   ns
   
}

##*****************************************************************************
##' Compute \code{NSGEV} quantiles.
##'
##' @title Compute \code{NSGEV} Quantiles
##'
##' @param x A \code{NSGEV} object. 
##'
##' @param probs Vector of probabilities.
##' 
##' @param data A data frame with the covariates. If \code{NULL}, the
##' data frame affached to the model given in \code{object} is used.
##'
##' @param psi Vecor of model coefficients.
##'
##' @param ... Not used yet.
##'
##' @return A matrix of quantiles, with one row for each observation
##' of \code{data} and one column by probability.
##'
##' @author Yves Deville
##'
##' @method quantile NSGEV
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##' ## model structure
##' ns <- NSGEV(formulas = list("loc" = ~ alpha + beta * t,
##'                             "scale" = ~ delta,
##'                             "shape" = ~ xi),
##'              data = df)
##' psi <- c("alpha" = 1, "beta" = 0.1, "delta" = 0.6, "xi" = 0.06)
##' q <- quantile(ns, psi = psi) 
##' matplot(x = df$t, y = q, type = "l", lty = 1:3, lwd = 2,
##'         main = "model quantiles")
##' legend("bottomright", legend = colnames(q), lty = 1:3, col = 1:3,
##'         lwd = 2)
quantile.NSGEV <- function(x, probs = c(0.90, 0.95, 0.99),
                           data = NULL,
                           psi = NULL, ...) {
    
    ## control (from  quantile.default)
    if (any(is.na(probs))) stop ("NA not allowed yet in 'probs'")
   
    eps <- 100 * .Machine$double.eps
    if (any(probs < eps | probs > 1 - eps))  stop("'probs' outside [0,1]")
    
    if (is.null(data)) data <- x$data
    if (is.null(psi)) psi <- x$estimate
    
    n <- nrow(data)
    theta <- psi2theta(model = x, psi = psi, data = data)
    quant <- array(NA, dim = c(n, length(probs)),
                   dimnames = list(rownames(data),
                       paste("Q", formatPerc(probs), sep = "")))

    for (i in seq_along(probs)) {
        quant[ , i] <- nieve::qGEV(probs[i], loc = theta[ , 1L],
                                   scale = theta[ , 2L],
                                   shape = theta[ , 3L])
    }
    attr(quant ,"p") <- probs
    
    quant
}

##*****************************************************************************
##' Compute NSGEV densities.
##'
##'
##' @title Compute NSGEV Densities
##'
##' @param x A \code{NSGEV} object.  Note that the name of this formal
##' is imposed by the S3 generic, it would otherwise probably have
##' been chosen as 'model'.
##'
##' @param xValue Vector of quantiles at which the GEV densities will be
##' evaluated. By default, a grid of value is found with coverage
##' probability \code{> 0.001} for each observation.
##' 
##' @param data A data frame with the covariates. If \code{NULL}, the
##' data frama affached to the model given in \code{x} is used.
##'
##' @param psi Vector of model coefficients.
##'
##' @param log Logical. If \code{TRUE} the log-density is returned.
##' 
##' @param ... Not used yet.
##'
##' @return A matrix of density values, with one row for each observation
##' of \code{data} and one column by quantile.
##'
##' @author Yves Deville
##'
##' @method density NSGEV
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##' ## model structure
##' ns <- NSGEV(formulas = list("loc" = ~ alpha + beta * t,
##'                             "scale" = ~ delta,
##'                             "shape" = ~ xi),
##'              data = df)
##' psi <- c("alpha" = 1, "beta" = 0.1, "delta" = 0.6, "xi" = 0.06)
##' d <- density(ns, psi = psi) 
##' matplot(x = attr(d, "x"), y = t(d), type = "l")
##' 
density.NSGEV <- function(x,
                          xValue = NULL,
                          data = NULL,
                          psi = NULL,
                          log = FALSE,
                          ...) {

    ## XXX change the 'quantile' method to optionnaly return more
    ## results to avoid recomputing them.
    if (is.null(xValue)) {
        q <- quantile(x, probs = c(5e-4, 1 - 5e-4),
                      data = data, psi = psi)
        r <- range(q)
        xValue <- seq(from = r[1], to = r[2], length.out = 200)
    } 
    
    if (is.null(data)) data <- x$data
    if (is.null(psi)) {
        psi <- x$estimate
    }
    
    n <- nrow(data)
    theta <- psi2theta(model = x, psi = psi, data = data)
    dens <- array(NA, dim = c(n, length(xValue)),
                  dimnames = list(rownames(data), NULL))
    
    for (i in seq_along(xValue)) {
        dens[ , i] <- nieve::dGEV(xValue[i], loc = theta[ , 1L],
                                  scale = theta[, 2L],
                                  shape = theta[, 3L], log = log)
    }
    attr(dens ,"x") <- xValue
    dens
}

##*****************************************************************************
##' Simulate paths from a NSGEV model conditional on the
##' covariates.
##'
##' @title Simulate Paths from a \code{NSGEV} Model Conditional on the
##' Covariates
##' 
##' @param object A NSGEV model.
##'
##' @param nsim Number of simulations (paths).
##'
##' @param seed Not used yet.
##'
##' @param data The data frame of covariates. By default, the model
##' data frame is used.
##'
##' @param psi Vector of model parameters. By default, the model
##' parameters are used.
##'
##' @param ... Not used yet.
##'
##' @return A matrix with \code{nrow(data)} rows and \code{nsim}
##' columns, each column representing a time.
##'
##' @author Yves Deville
##'
##' @seealso \code{\link[nieve]{rGEV}} to simulate from varying GEV
##' parameters.
##'
##' @method simulate NSGEV
##' @export
##' 
##' @examples
##' df <- data.frame(t = 1:10)
##' ## model structure
##' ns <- NSGEV(formulas = list("loc" = ~ alpha + beta * t,
##'                             "scale" = ~ delta,
##'                             "shape" = ~ xi),
##'              data = df)
##' df.new <- data.frame(t = 11:20)
##' psi <- c("alpha" = 1, "beta" = 0.1, "delta" = 0.6, "xi" = 0.06)
##' ysim <- simulate(ns, nsim = 20, psi = psi)
##' ysim.new <- simulate(ns, nsim = 20, data = df.new, psi = psi)
##' matplot(df$t, ysim, type = "l", xlab = "t",
##'         col = "orange", xlim = c(1, 20))
##' matlines(df.new$t, ysim.new, type = "l", col = "firebrick")
simulate.NSGEV <- function(object, nsim = 1, seed = NULL,
                           data = NULL, psi = NULL, ...) {
   
   if (is.null(data)) data <- object$data
   if (is.null(psi)) psi <- object$estimate
   
   n <- nrow(data)
   theta <- psi2theta(model = object, psi = psi, data = data) 
    sim <- nieve::rGEV(nsim, loc = theta[ , 1L],
                       scale = theta[ , 2L],
                       shape = theta[ , 3L])
   dimnames(sim) <- list(rownames(data), paste("sim", 1L:nsim, sep = ""))
   
   sim
   
}




