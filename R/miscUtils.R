## ***********************************************************************
##' Translate a standard formula into a NSGEV formula.
##'
##' @title Translate a Formula into a NSGEV Formula
##' 
##' @param formula The formula to be translated.
##'
##' @param parnm The name of the parameter to use, i.e.
##' "mu,""sigma" or "xi".
##'
##' @return A formula with the parameters displayed, suitably
##' named. The parameter names are attached to the result as an
##' attribute.
##'
##' @author Yves Deville
##'
##' @export
##' 
##' @examples
##' transFormula(~time, parnm = "mu") 
##' 
transFormula <- function(formula, parnm) {
    vnms <- all.vars(formula)
    nv <- length(vnms)
    if (nv) {
        parnms <- paste(parnm, vnms, sep = "_")
        text <- paste(paste(parnms, rep("*", nv), vnms), collapse = " + ")
    } else {
        text <- parnms <- character(0)
    }
    parnm0 <- paste(parnm, "0", sep = "_")
    parnms <- c(parnm0, parnms)
    text <- paste(c(parnm0, text), collapse = " + ")
    res <- as.formula(sprintf(" ~ %s", text))
    attr(res, "parNames") <- parnms
    res
}

## ***********************************************************************
##' Generates random names.
##'
##' This function is a utility function for tests concerning the
##' parsing of formulas.
##' 
##' @title Generates Random Names
##'
##' @param n Number of names.
##'
##' @param nchar Maximal number of characters for each name.
##' 
##' @return A character vector of valid names.
##'
##' @export
##'  
##' @examples
##' set.seed(31415)
##' rNames(4)
##' rNames(4, nchar = 6)
rNames <- function(n, nchar = 3L) {
    df <- expand.grid(LETTERS, LETTERS, LETTERS)
    LET <- apply(X = df[sample(1:676, n), ],
                 MARGIN = 1L,
                 FUN = function(x) x[sample(1L:3L,
                                            size = sample(1L:3L, size = 1))])
    paste0(sapply(LET, paste0, collapse = ""), "")
}

## copy of stats:::format_perc

formatPerc <- function (x,
                         digits = max(2L, getOption("digits")),
                         probability = TRUE, 
                         use.fC = length(x) < 100,
                         ...) {
    if (length(x)) {
        if (probability)  x <- 100 * x
        paste0(if (use.fC) 
                   formatC(x, format = "fg", width = 1, digits = digits)
               else format(x, trim = TRUE, digits = digits, ...), "%")
    }
    else character(0)
}

## pos <- function(date, br = "1970-01-01") {
##     res <- as.numeric(as.Date(date) - as.Date(br))
##     res[res < 0] <- 0
##     res
## }

## ***********************************************************************
##' Find block duration in years.
##'
##' This is simply a wrapper for the \code{diff} method.
##' 
##' @title Find Block Duration 
##' 
##' @param date A vector with class \code{"Date"} or that can be
##' coerced to this class.
##'
##' @return A duration in years.
##'
##' @note For most applications of Time Varying GEV models with
##' \code{TVGEV}, the block duration will be one year. Yet it can be
##' sometimes needed to form blocks with longer duration (e.g. two
##' years, five years).
##'
##' @export
##' 
##' @examples
##'
##' date <- as.Date(sprintf("%4d-01-01", TXMax_Dijon$Year))
##' n <- length(date)
##' blockDuration(date)
##' blockDuration(date[seq(from = 1, to = n, by = 2)])
blockDuration <- function(date) {

    date <- as.Date(date)
    d <- as.numeric(diff(date), units = "days")
    
    if (sd(d) / mean(d) > 0.005) {
        warning("unevenly spaced dates")
    }

    round(mean(d) / 365.25, digits = 3)
    
}


formatLevel <- function(level) {
    paste(gsub("\\.0$", "", sprintf("%4.1f", 100 * level)), "%", sep = "")
}

## copy of stats:::format_perc

formatPerc <- function(x,
                       digits = max(2L, getOption("digits")),
                       probability = TRUE, 
                       use.fC = length(x) < 100,
                       ...) {
    if (length(x)) {
        if (probability)  x <- 100 * x
        paste0(if (use.fC) 
                   formatC(x, format = "fg", width = 1, digits = digits)
               else format(x, trim = TRUE, digits = digits, ...), "%")
    }
    else character(0)
}

## pos <- function(date, br = "1970-01-01") {
##     res <- as.numeric(as.Date(date) - as.Date(br))
##     res[res < 0] <- 0
##     res
## }

# ***********************************************************************
##' Select (nearly) a given number of dates.
##' 
##' @title Select (Nearly) a Given Number of Dates
##'
##' @param date A vector with class\code{"Date"} or that can be coerced
##' to this class.
##'
##' @param n Target number of selected dates.
##'
##' @return A vector of 'round' dates having nearly length \code{n}
##' and covering approximately the range of \code{date}.
##'
##' @export
##' 
selectDate <- function(date, n = 3L) {

    date <- as.Date(date)
    
    rd <- range(date)
    W <- as.numeric(diff(rd, unit = "days") ) / 365
    ints <- c(1, 5, 10, 20, 25, 50, 100, 200, 500)
    nw <- round(W / ints)
    n - 1
    i <- which.min(abs(outer(n - 1, nw, "-")))
    w <- ints[i]
    y <- as.numeric(format(rd, "%Y"))
    
    if (W > 30) {
        y1 <- 10 * floor(y[1] / 10)
        ys <- seq(from = y1, by = w, length.out = n)
    }

    as.Date(sprintf("%4d-01-01", ys))
}

## *****************************************************************************
##' Choose an opacity level.
##'
##' @title Choose an Opacity Level
##'
##' @param n Vector of integers.
##' 
##' @return Vector of values for the opacity "alpha". Value 0 is for a
##' fully transparent colour, value 1 is for a fully opaque colour.
##'
##' @export
##' 
opacity <- function(n) {
    exp(- log(n)^2 / 16)
}

## *****************************************************************************
##' Distance between two lines.
##' 
##' This function is used to build diagnostics during the
##' construction of profile-likelihood confidence intervals. So the
##' function is not intended to be used by itself.
##' 
##' @title Distance Between Two Lines
##' 
##' @param x1 A numeric vector with positive norm.
##' 
##' @param x2 A numeric vector with positive norm and with the same
##' length as \code{x1}.
##' 
##' @return The distance.
##'
##' @details The vectors \code{x1} and \code{x2} are scaled to have an
##' unit Euclidean norm and to have the same direction. The distance
##' between the lines is simply the distance between the unit vectors
##' on the unit sphere.
##'
##' @note When a constrained optimisation is used for the
##' determination of the lower or upper bound of the interval, we have
##' to check that the gradient of the objective and that of the
##' constraint are colinear, which is performed with this function.
##'
##' @export
##' 
distLines <- function(x1, x2) {
    if (length(x1) != length(x2)) {
        stop("'x1' and 'x2' must have the same length")
    }
    if (any(is.na(x1)) || any(is.na(x2))) {
        return(NA)
    }
    n1 <- sqrt(sum(x1^2))
    n2 <- sqrt(sum(x2^2))
    if ((n1 <= 1e-6) || (n2 <= 1e-6)) {
        ## warning("'x1' and 'x2' should have Euclidean norm > 1e-6")
        return(NA)
    }
    x1 <- x1 / n1
    x2 <- x2 / n2
    i <- which.max(abs(x1))
    if (x1[i] * x2[i] < 0.0) x2 <- -x2
    s <- sum(x1 * x2)
    if (s > 1 - 1e-9) return(0.0)
    acos(s)
}


    
