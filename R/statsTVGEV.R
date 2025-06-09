
## *****************************************************************************
##' Compute quantiles for a \code{TVGEV} object representing a
##' time-varying model.
##'
##' @title  Compute Quantiles for a \code{TVGEV} Object
##'
##' @param x A \code{TVGEV} object representing a time-varying GEV
##' model.
##'
##' @param probs Vector of probabilities.
##'
##' @param date An object that can be coerced to the class
##' \code{"Date"}. The default \code{NULL} corresponds to using
##' the dates attached to \code{object}.
##'
##' @param psi An optional vector of parameters for \code{object}.
##' 
##' @param ... Not used yet.
##'
##' @return An object with class \code{quantile.TVGEV} inheriting from
##' \code{"matrix"}. This is essentially a matrix with one column by
##' quantile, but it has a \code{"date"} attributes than can be used
##' for plotting.
##'
##' @importFrom stats quantile
##' @method quantile TVGEV
##' @export
##' 
##' @examples
##' example(TVGEV)
##' q <- quantile(res1)
##' autoplot(q)
##' 
quantile.TVGEV <- function(x, probs = c(0.90, 0.95, 0.99),
                           date = NULL,
                           psi = NULL, ...) {

    if (length(x$TSVars)) {
        stop("'x' includes TSVars. It can not be used for now")
    }
    
    ## control (from  quantile.default)
    if (any(is.na(probs))) stop ("NA not allowed yet in 'probs'")
    probs <- sort(probs)
    
    eps <- 100 * .Machine$double.eps
    if (any(probs < -eps | probs > 1 + eps))  stop("'probs' outside [0,1]")

    if (is.null(psi)) psi <- x$estimate

    newdate <- checkNewdata.TVGEV(object = x, newdata = date)
    fDate <- format(newdate[[x$date]])
    theta <- psi2theta(model = x, psi = psi, date = newdate)
    n <- nrow(theta)
    fProbs <- formatPerc(probs)
    
    quant <- array(NA, dim = c(n, length(probs)),
                   dimnames = list(date = fDate,
                       paste("Q", fProbs, sep = "")))

    for (i in seq_along(probs)) {
        quant[ , i] <- nieve::qGEV(p = probs[i], loc = theta[ , 1L],
                                   scale = theta[ , 2L],
                                   shape = theta[ , 3L])
    }
    
    attr(quant, "date") <- fDate
    attr(quant, "collabels") <- fProbs
    attr(quant, "label") <- "quantile"
    class(quant) <- c("bts", "matrix")
    return(quant)
    
}
## *************************************************************************
##' Compute TVGEV densities or Cumulative Distribution Functions.
##'
##' @rdname density.TVGEV
##'
##' @aliases density.TVGEV
##' @aliases cdf.TVGEV
##' 
##' @title Compute TVGEV Densities or Cumulative Distribution
##' Functions
##'
##' @param x A \code{TVGEV} object.
##'
##' @param xValue Vector of quantiles at which the GEV densities will
##' be evaluated. By default, a grid of value is found with coverage
##' probability > 0.001 for each observation.
##' 
##' @param date An object that can be coerced to the
##' class\code{"Date"} giving the date of the blocks for which the
##' density will be evaluated. By default, three dates are selected
##' in the date vector attached to \code{x}.
##'
##' @param psi Vector of model coefficients. By default, the vector of
##' estimated coefficients in \code{x} is used.
##'
##' @param log Logical. If TRUE the log-density is returned.
##'
##' @param ... Other arguments to be passed to
##'     \code{\link[nieve]{dGEV}} or \code{\link[nieve]{pGEV}}.
##'
##' @return An oject with class \code{"bfts"}. This is mainly a matrix
##' with one row by date. Rather than providing a time-plot of each
##' column (as would be done for a \code{"bts"} object), the
##' \code{plot} method plots the density or cdf function for a small
##' number of dates. Unless the number of dates is \code{1}, the
##' functions are plotted with the x-y axes flipped in order to enhance
##' the time-varying feature of the model.
##' 
##' @seealso \code{\link{GEV}} for the density and cdf of the GEV
##' distribution, \code{\link{plot.predict.TVGEV}} for the Return Level
##' plot.
##'
##' @importFrom stats density
##' @method density TVGEV
##' @export
##' 
##' @examples
##' example(TVGEV)
##' d <- density(res1, date = c("2000-01-01", "2020-01-01", "2040-01-01"))
##' F <- cdf(res1, date = c("2000-01-01", "2020-01-01", "2040-01-01"))
##' plot(d, fill = TRUE)
##' plot(F)
##' 
density.TVGEV <- function(x, xValue = NULL,
                          date = NULL,
                          psi = NULL,
                          log = FALSE, ...) {

    if (length(x$TSVars)) {
        stop("'x' includes TSVars. It can not be used for now")
    }

    ## avoid using two many dates
    if (is.null(date)) date <- selectDate(x$fDate)
    newdate <- checkNewdata.TVGEV(object = x, newdata = date)
    fDate <- format(newdate[[x$date]])
    n <- nrow(newdate)
    
    if (is.null(psi)) psi <- x$estimate
    theta <- psi2theta(model = x, psi = psi, date = newdate)

    if (is.null(xValue)) {
        q <- quantile(x, probs = c(0.001, 0.999), date = newdate)
        rq <- range(q)
        xValue <- seq(from = rq[1L], to = rq[2L], length.out = 200)
    }
    fxValue <- format(xValue)
    dens <- array(NA, dim = c(n, length(xValue)),
                  dimnames = list(date = fDate,
                                  paste("d", fxValue, sep = "")))
    
    for (i in seq_along(xValue)) {
        dens[ , i] <- nieve::dGEV(xValue[i], loc = theta[ , 1L],
                                  scale = theta[, 2L],
                                  shape = theta[, 3L], ...)
    }
    
    attr(dens, "x") <- xValue
    rownames(dens) <- attr(dens, "date") <- fDate
    attr(dens, "collabels") <- fxValue
    attr(dens, "label") <- "density"
    class(dens) <- c("bfts", "bts", "matrix")
    return(dens)


}

## *************************************************************************
##' @rdname density.TVGEV
##' 
##' @param qValue xValue Vector of probabilities at which the GEV
##' densities will be evaluated. By default, a grid of value from
##' \code{0.0} to \code{1.0} is used. for each observation.
##'
##' @method cdf TVGEV
##' @export
##' 
cdf.TVGEV <- function(x,
                      qValue = NULL,
                      date = NULL,
                      psi = NULL,
                      log = FALSE, ...) {

    if (length(x$TSVars)) {
        stop("'object' includes TSVars. It can not be used for now")
    }

    ## avoid using two many dates
    if (is.null(date)) date <- selectDate(x$fDate)
    newdate <- checkNewdata.TVGEV(object = x, newdata = date)
    fDate <- format(newdate[[x$date]])
    n <- nrow(newdate)
    
    if (is.null(psi)) psi <- x$estimate
  
    theta <- psi2theta(model = x, psi = psi, date = newdate)
    
    if (is.null(qValue)) {
        q <- quantile(x, probs = c(0.005, 0.995), date = newdate)
        rq <- range(q)
        qValue <- seq(from = rq[1L], to = rq[2L], length.out = 200)
    }
    
    fqValue <- format(qValue)
    cdf <- array(NA, dim = c(n, length(qValue)),
                 dimnames = list(date = format(date),
                     paste("d", fqValue, sep = "")))
    
    for (i in seq_along(qValue)) {
        cdf[ , i] <- nieve::pGEV(qValue[i], loc = theta[ , 1L],
                                 scale = theta[, 2L],
                                 shape = theta[, 3L], ...)
    }
    
    attr(cdf, "x") <- qValue
    rownames(cdf) <- attr(cdf, "date") <- fDate
    attr(cdf, "collabels") <- fqValue
    attr(cdf, "label") <- "cdf"
    class(cdf) <- c("bfts", "bts", "matrix")
    return(cdf)

}

## *************************************************************************
##' Expectation or moments of a \code{TVGEV} model.
##'
##' @name mean.TVGEV
##' @rdname mean.TVGEV
##' @aliases mean.TVGEV
##' @aliases moment.TVGEV
##'
##' @title Expectation or Moments of a \code{TVGEV} Model Object
##'
##' @param x An \code{TVGEV} object.
##'
##' @param date An object that can be coerced to the \code{"Date"}
##' class. By default the date in \code{x} is used. 
##'
##' @param psi Vector of parameters. By default the vector of the
##' parameters for the model \code{x} is used.
##'
##' @param ... Not used yet.
##'
##' @return An object with class \code{"bts"} representing a block
##' time series with its value set to the expectation or the moment of
##' the GEV distribution.
##'
##' @seealso \code{\link{quantile.TVGEV}}
##' to compute quantiles of the time-varying distribution.
##'
##' @method mean TVGEV
##' @export
##' 
##' @examples
##' example(TVGEV)
##' autoplot(mean(res1))
mean.TVGEV <- function(x,
                       date = NULL,
                       psi = NULL, ...) {

    if (length(x$TSVars)) {
        stop("'x' includes TSVars. It can not be used for now")
    }
     
    ## Euler-Mascheroni constant
    gam <- - digamma(1) 
    if (is.null(psi)) psi <- x$estimate

    newdate <- checkNewdata.TVGEV(object = x, newdata = date)
    fDate <- format(newdate[[x$date]])
    n <- nrow(newdate)
   
    theta <- psi2theta(model = x, psi = psi, date = newdate)

    E <-  theta[ , 1, drop = FALSE] + gam * theta[ , 2, drop = FALSE]
    colnames(E) <- "expectation"
    ind <- theta[ , 3] != 0.0
    E[ind, 1L] <- theta[ , 1] + theta[ind, 2] * (gamma(1 - theta[ind, 3]) - 1.0 ) / theta[ind, 3]
    
    rownames(E) <- attr(E, "date") <- fDate
    attr(E, "collabels") <- "expectation"
    attr(E, "label") <- "expectation"
    class(E) <- c("bts", "matrix")
    return(E)
}


## *************************************************************************
##' @rdname mean.TVGEV
##' @name moment.TVGEV
##'
##' @param which Description of the moment. For now only the value
##' \code{"variance"} is accepted, but a description for the 3-rd and 4-th
##' order moments should be accepted in the future.
##'
##' @method moment TVGEV
##' @export 
moment.TVGEV <- function(x,
                         which = "variance",
                         date = NULL,
                         psi = NULL, ...) {
    if (length(x$TSVars)) {
        stop("'x' includes TSVars. It can not be used for now")
    }
 
    if (which != "variance") {
        stop("'which' can only be \"variance\" for now.")
    }
    
    ## Euler-Mascheroni constant
    gam <- - digamma(1)
    
    if (is.null(psi)) psi <- x$estimate
    
    newdate <- checkNewdata.TVGEV(object = x, newdata = date)
    fDate <- format(newdate[[x$date]])
    n <- nrow(newdate)
  
    theta <- psi2theta(model = x, psi = psi, date = newdate)
    
    M <-  array(pi * pi / 6, dim = c(n, 1L),
                dimnames = list(date = fDate,
                    "moment" = "variance"))
    ## \xi >= 1/2
    ind <- (theta[ , 3L] >= 0.5)
    M[ind, ] <- Inf

    ## \xi = 0
    ind <- theta[ , 3L] != 0.0
    M[ind, 1L] <- theta[ind, 2L]^2 *
        (gamma(1 - 2 * theta[ind, 3L]) - gamma(1 - theta[ind, 3L])^2 ) /
            theta[ind, 3] / theta[ind, 3]
    
    attr(M, "date") <- fDate
    attr(M, "collabels") <- "variance"
    attr(M, "label") <- "variance"
    class(M) <- c("bts", "matrix")
    return(M)
}


