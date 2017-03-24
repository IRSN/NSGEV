##' Build a design matrix to describe breaks in a linear trend.
##'
##' The design matrix is a Truncated Power Basis of splines. It
##' can be extended to a degree differing from the default value \eqn{1}.
##' 
##' @title Design Matrix with Breaks in a Linear Trend
##' 
##' @param date A vector with class \code{"Date"} or a vector that can
##' be coerced to the \code{"Date"} class, typically an unambigous
##' character vector.
##'
##' @param degree The degree used for spline functions.
##'
##' @param origin Not used yet.
##'
##' @param breaks A vector describing the breaks. It will be coerced
##' to the \code{"Date"} class, and can be e.g. an unambigous
##' character vector.
##'
##' @param constant Logical. If \code{TRUE} a column with constant
##' unit value is added.
##'
##' @return A matrix with the value of the basis functions as its
##' columns.
##'
##' @details
##'
##' Consider the date as a numeric variable \eqn{t}, let \eqn{d} be
##' the degree.  The basis consist in
##'
##' \enumerate{
##'
##' \item the powers \eqn{t^0 = 1}, \eqn{t^1}, \eqn{\dots}{...},
##' \eqn{t^d} where \eqn{d} is the given degree.
##'
##' \item the truncated powers \eqn{(t - \tau_k)^d_+} where \eqn{\tau_k}
##' is the break with number \eqn{k} and \eqn{x_+} is the positive
##' part of \eqn{x} (so \eqn{x_+ = 0} when \eqn{x < 0}).
##'
##' }
##'
##' The power zero will be discarded when \code{contant} is
##' \code{FALSE}.
##'
##' @note When the argument \code{breaks} has length zero, the
##' function returns a polynomial basis with the constant possibly
##' removed.
##' 
##' @author Yves Deville
##'
##' @examples
##' 
##' date <- as.Date(sprintf("%4d-01-01", 1921:2020)) 
##' X1 <- breaksX(date = date, breaks = c("1970-01-01", "1990-01-01"))
##' plot(date, X1[ , 1], ylim = range(X1), type = "n", ylab = "")
##' matlines(date, X1, type = "l") 
breaksX <- function(date,
                    degree = 1L,
                    origin = NULL,
                    breaks = NULL,
                    constant = FALSE) {
    
    breaks <- as.Date(breaks)
    tBreaks <- as.numeric(breaks)  / 365.25
    t <- as.numeric(date) / 365.25

    if (constant) { 
        X <- outer(t, 0:degree, FUN = "^")
        colnames(X) <- c("cst", paste("t", 1:degree, sep = ""))
    } else {
        X <- outer(t, 1:degree, FUN = "^")
        colnames(X) <- paste("t", 1:degree, sep = "")
    }
    
    if (!is.null(breaks)) {
        dif <- outer(X = t, Y = tBreaks, FUN = "-")
        colnames(dif) <- paste(paste("t", degree, sep = ""),
                               format(breaks, "%Y"), sep = "_")
        dif[dif < 0] <- 0
        if (degree != 1) {
            dif[dif > 0] <- dif[dif > 0]^degree
        }
        X <- cbind(X, dif)
    }

    X
    
}
