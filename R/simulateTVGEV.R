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
##' @importFrom stats simulate
##' @method simulate TVGEV
##' @export
##'  
##' @examples
##' example(TVGEV)
##' sim <- simulate(res2, nsim = 200)
##' plot(sim)
##' 
simulate.TVGEV <- function (object, nsim = 1, seed = NULL,
                            newdate = NULL, psi = NULL,  ...) {
    
    theta <- psi2theta(model = object, psi = psi, date = newdate)
    
    sim <- nieve::rGEV(nsim, loc = theta[ , 1L],
                       scale = theta[ , 2L], 
                       shape = theta[ , 3L])
    
    dimnames(sim) <- list(rownames(theta),
                          paste("sim", 1L:nsim,  sep = ""))

    if (is.null(newdate)) newdate <- object$data[ , object$date]
    
    attr(sim, "date") <- newdate
    
    ## class(sim) <- "simulate.TVGEV"
    class(sim) <- c("bts", "matrix")
    sim

}

## *************************************************************************
## THIS METHOD NO LONGER EXISTS SINCE THE RESULT OF THE 'simulate' METHOD
## NOW HAS CLASS "bts"

##' 
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
##' method plot simulate.TVGEV
##' export
##' @noRd
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

