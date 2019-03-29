## ****************************************************************************
##' Build a design matrix to describe breaks in a linear trend.
##'
##' The design matrix is a Truncated Power Basis of splines. It can be
##' extended to a degree differing from the default value \eqn{1}.
##' The provided date is converted into a numeric variable \eqn{t}
##' corresponding to the number of years from the origin. The basis
##' functions are functions of this variable \eqn{t}. A year is
##' defined as \code{365.25} days, so \eqn{t} can not take only
##' integer values.
##' 
##' @title Design Matrix with Breaks in a Linear Trend
##' 
##' @param date A vector with class \code{"Date"} or a vector that can
##' be coerced to the \code{"Date"} class, typically an unambiguous
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
##' @return An object of class \code{"bts"} inheriting from
##' \code{"matrix"}. This is essentially a numeric matrix with the
##' value of the basis functions as its columns. The rows are in
##' correspondence with the elements of the \code{date} vector given
##' on input.
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
##' The power zero will be discarded when \code{constant} is
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
##' plot(X1) 
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

    attr(X, "date") <- rownames(X) <- format(date, "%Y-%m-%d")
    class(X) <- c("bts", "matrix")
    
    X
    
}

## ****************************************************************************
##' Build a design matrix for polynomial regression.
##'
##' The provided date is converted into a numeric variable \eqn{t}
##' corresponding to the number of years from the origin. The basis
##' functions are functions of this variable \eqn{t}. A year is
##' defined as \code{365.25} days, so \eqn{t} can not take only
##' integer values.
##' 
##' @title Design Matrix for Polynomial Regression
##' 
##' @param date A vector with class \code{"Date"} or that can be
##' coerced to this class.
##' 
##' @param degree The maximal degree \eqn{d}. The basis contains
##' \eqn{d + 1} functions \eqn{t^0}, \eqn{t^1}, \dots, \eqn{t^d}.
##'
##' @param origin Optional vector of length 1 with class \code{"Date"}
##' or character that can be coerced to this class. Gives the origin
##' for the transformation of dates into numeric values. The default
##' value is \code{1970-01-01}, see \code{\link{Date}}.
##'
##' @return An object with class \code{"bts"} inheriting from
##' \code{"matrix"}. This is essentially a numeric matrix having the
##' value of the basis functions as its columns. Each row corresponds
##' to an element of the given \code{date} vector. By convention, the
##' columns are named \code{Cst} then \code{t1}, \code{t2}, \dots.
##' 
##' @note This function is intended to provide regressors for Time
##' Varying GEV models corresponding to \code{TVGEV} objects. Since
##' the response then usually corresponds to block maxima, the date
##' variable gives either the begining or the end of each block.
##'
##' @seealso \code{breaksX} for a basis of Truncated Power functions.
##' and \code{TVGEV} to for Time Varying GEV models.
##'
##' @section Caution: As a general rule, it is a good practice to
##' choose \code{origin} close to the centre of the time period. This
##' is especially important when a polynomial of degree \eqn{>= 2} is
##' used since problems of convergence are likely to occur otherwise.
##' 
##' @examples
##'
##' date <- seq(from = as.Date("1996-01-01"),
##'             to = as.Date("2016-12-31"), by = "years")
##'
##' X <- polynomX(date)
##' plot(X)
polynomX <- function(date, degree = 2, origin = NULL) {

    if (degree < 1) stop("'degree' must be >= 1")
    ## dtNum <- as.numeric(diff(head(date)), units = "days")
    
    mc <- match.call()
    date <- as.Date(date)
    n <- length(date)

    if (!is.null(origin)) {
        origin <- as.Date(origin)
        tt <- as.numeric(date - origin) / 365.25
    } else {
        tt <- as.numeric(date) / 365.25
    }
    X <- outer(tt, 0:degree, FUN = "^")
    colnames(X) <- c("Cst", paste("t", 1:degree, sep = ""))
    rownames(X) <- attr(X, "date") <- format(date, "%Y-%m-%d")
    attr(X, "origin") <- origin
    class(X) <- c("bts", "matrix")
    X
}

## ****************************************************************************
##' Build a design matrix for a natural spline basis.
##'
##' The design matrix contains a basis to represent a natural spline
##' with possibly chosen knots if needed. By default the boundary
##' knots are chosen as the range of the numeric translation of
##' \code{date}. The basis generates a space of cubic splines with
##' their second order derivative vanishing at the two boundary
##' knots. Therefore, the splines behave as a polynomials of degree
##' one outside of the boundary knots.
##' 
##' @title Design Matrix for a Natural Spline Basis
##' 
##' @param date A vector with class \code{"Date"} or a vector that can
##' be coerced to the \code{"Date"} class, typically an unambiguous
##' character vector.
##'
##' @param origin Optional vector of length 1 with class \code{"Date"}
##' or character that can be coerced to this class. Gives the origin
##' for the transformation of dates into numeric values. The default
##' value is \code{1970-01-01}, see \code{\link{Date}}.
##'
##' @param knots A vector with class \code{"Date"} or a vector that
##' can be coerced to the \code{"Date"} class, giving the location of
##' the knots as in \code{\link[splines]{ns}}.
##' 
##' @param boundaryKnots A vector with class \code{"Date"} or a
##' vector that can be coerced to the \code{"Date"} class, giving the
##' location of the boundary knots in \code{\link[splines]{ns}}.
##'
##' @param constant Logical passed as \code{intercept} to
##' \code{\link[splines]{ns}}.
##' 
##' @return An object with class \code{"bts"} inheriting from
##' \code{"matrix"}. This is essentially a numeric matrix having the
##' value of the basis functions as its columns. Each row corresponds
##' to an element of the given \code{date} vector. By convention, the
##' columns are named \code{"ns1"}, \code{"ns2"}, \dots The number of
##' colums is equal to \code{df} is provided or if knots was supplied,
##' to \code{length(knots) + 1 + constant}.
##'
##' @seealso The \code{\link[splines]{ns}} function used to compute
##' the basis, and \code{\link{breaksX}}, \code{\link{polynomX}} for
##' alternative bases.
##'
##' @examples
##' date <- as.Date(sprintf("%4d-01-01", 1921:2020))
##' X1 <- natSplineX(date = date, knots = "1950-01-01",
##'                  boundaryKnots = c("1920-01-01", "2021-01-01"))
##' plot(X1)
##' 
natSplineX <- function(date,
                       origin = NULL,
                       knots,
                       boundaryKnots,
                       constant = TRUE) {

    t <- as.numeric(date) / 365.25
    
    knots <- as.Date(knots)
    tKnots <- as.numeric(knots) / 365.25
    
    boundaryKnots <- as.Date(boundaryKnots)
    tboundaryKnots <- as.numeric(boundaryKnots)  / 365.25
    
    X <- splines::ns(t, df = df,  knots = tKnots, intercept = constant,
            Boundary.knots = tboundaryKnots)
    
    colnames(X) <- paste("ns", 1:ncol(X), sep = "")
    rownames(X) <- attr(X, "date") <- format(date, "%Y-%m-%d")

    class(X) <- c("bts", "matrix")
    
    X
    
}

