
##' Coerce into a \code{NSGEV} object.
##'
##' @title Coerce into a \code{NSGEV} Object
##'
##' @param x Object to be coerced.
##'
##' @param ... Special arguments that may be required for some
##' classes.
##' 
##' @return An object with (S3) class \code{"NSGEV"} representing
##' a Non-Stationary GEV model.
##' 
as.NSGEV <- function(x, ...) {
    UseMethod("as.NSGEV")
}

##' Coerce a 'gev.fit' object into a 'NSGEV' object.
##'
##' @title Coerce a 'gev.fit' Object into a 'NSGEV' Object
##' 
##' @param x The object to be coerced.
##'
##' @param data The matrix of covariates that was used.
##' 
##' @param ... Not used yet.
##'
##' @return An object with class \code{"NSGEV"}.
##'
##' @note A \code{gev.fit} object does not embed the set of covariates
##' used nor a even a call. Thus the set of covariates must be given.
##' To be compliant with \code{NSGEV}, \emph{this data set must be a data frame
##' or a matrix with colnames}.
##' 
##' @author Yves Deville
##'
##' @examples
##' require(ismev)
##' ## number of observations and covariates
##' n <- 20; m <- 4 
##' 
##' ## generate a matrix covariates
##' set.seed(1234)
##' dat <- matrix(runif(n * m), nrow = n)
##' colnames(dat) <- rNames(m)
##'
##' ## response
##' y <- drop(rGEV(n))
##'
##' ## fit 
##' fit <- gev.fit(xdat = y,             ## response     
##'                ydat = dat,           ## matrix of covariates
##'                mul = 1:2, sigl = 4,  ## indices of covariates
##'                shl = NULL,           
##'                mulink = exp)         ## optional inverse-link
##'
##' ## coerce
##' ns <- as.NSGEV(fit, data = dat)
##' 
as.NSGEV.gev.fit <- function(x, data, ...) {

    data <- as.data.frame(data)
    
    parNames.GEV <- c("loc", "scale", "shape")
    symNames.GEV <- c("mu", "sigma", "xi")
    
    forms <- list()
    parNames <- character(0)
    for (i in 1L:3L) {
        if (is.null(x$model[[i]])) {
            forms[[parNames.GEV[i]]] <- symNames.GEV[i]
            parNames <- c(parNames, symNames.GEV[i])
        } else {
            parNames <- c(parNames, symNames.GEV[i])
            nm <- colnames(data)[x$model[[i]]]
            nterms <- length(nm)
            pars <- paste(symNames.GEV[i], nm, sep = "_")
            parNames <- c(parNames, pars)
            text <- paste(pars, rep("*", nterms), nm, sep = "")
            text <- paste(text, collapse = " + ")
            text <- paste(symNames.GEV[i], text, sep = " + ")
            forms[[parNames.GEV[i]]] <- text
        }
    }
    
    nc <- nchar(x$link)
    links <- unlist(strsplit(substr(x$link, start = 3, stop = nc - 1), split = ","))
    links <- gsub("^\\s+", "", links)
    
    for (i in 1:3) {
        if (links[i] != "identity") {
            forms[[i]] <- as.formula(sprintf("~ %s(%s)", links[i], forms[[i]]))
        } else {
            forms[[i]] <- as.formula(sprintf(" ~ %s", forms[[i]]))
        }
    }
    
    psi <- x$mle
    names(psi) <- parNames
    
    ns <- NSGEV(formulas = forms, data = data, psi = psi)

}

##' Coerce a 'fevd' object into a 'NSGEV' object.
##'
##' @title Coerce a 'fevd' Object into a 'NSGEV' object
##'
##' @param x The object to be coerced. Must have class \code{"fevd"} and
##' be of type \code{"GEV"}.
##'
##' @param ... Not used yet.
##'
##' @return An object with class \code{"NSGEV"}.
##'
##' @author Yves Deville
##'
##' @section Caution: For now, only simple formulas can be used. 
##' 
##' @examples
##' require(extRemes)
##' 
##' ## see the examples for extRemes::fevd.
##' data(PORTw)
##'
##' ## fit a GEV model
##' fit <- fevd(x = TMX1, data = PORTw,
##'             location.fun = ~AOindex, scale.fun = ~AOindex,
##'             units = "deg C")
##'
##' ## coerce
##' ns <- as.NSGEV(fit)
##' 
as.NSGEV.fevd <- function(x, ...) {

    ## check parnames with x$parnames and rebuilt
    parNames.GEV <- c("loc", "scale", "shape")
    symNames.GEV <- c("mu", "sigma", "xi")
    longParNames.GEV <- c("location", "scale", "shape")

    if (x$type != "GEV") stop("'x' must be of type \"GEV\"")

    forms <- list()
    parNames <- character(0)
    for (i in  seq_along(longParNames.GEV)) {
        nma <-  sprintf("%s.fun", longParNames.GEV[[i]])
        tf <- transFormula(formula = x$call[[nma]],
                         parnm = symNames.GEV[[i]])
        forms[[parNames.GEV[[i]]]] <- tf
        parNames <- c(parNames, attr(tf, "parNames"))
    }
    psi <- x$results$par
    names(psi) <- parNames
        
    ns <- NSGEV(formulas = forms,
                data = x$cov.data,
                ## response = x$x,
                psi = psi,
                est = "none")
    ns$response <- x$x
    ns$negLogLik <- x$results$value
    vcov <- parcov.fevd(x)
    ## patch because some problems are met with hessian evaluation is
    ## package 'extRemes'. Use 'numDeriv' instead.
    if (is.null(vcov)) {
        hessian <- hessian(negLogLik, psi, model = ns)
        vcov <- solve(hessian)
    }
    rownames(vcov) <- colnames(vcov) <- parNames
    ns$vcov <- vcov
    ns

}
