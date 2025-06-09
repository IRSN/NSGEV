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
##' This function is used to select the dates for a conditional
##' prediction it those are not given by the user.
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

## 
## Copied from Renext
##
##' @importFrom grDevices col2rgb
##' 
translude <- function (colors, alpha = 0.6) {
    L <- pmax(length(colors), length(alpha))
    colors <- rep(colors, length.out = L)
    alpha <- rep(alpha, length.out = L)
    rgb <- as.matrix(col2rgb(colors) / 255)
    colors2 <- rgb(red = rgb["red", ], green = rgb["green", ], 
                   blue = rgb["blue", ], alpha = alpha)
}

##' Test that a provided object containing derivatives is in good
##' accordance with another similar object containing numeric
##' derivatives.
##'
##' One could think of requiring a small absolute difference or a
##' small relative difference between the results. However, in
##' practice it seems that the absolute difference can be large even
##' with correct derivatives because these contain large values. Also
##' the relative difference can be large if the derivatives contain
##' values that are small in absolute value. However, if one of these
##' two condition is fulfilled, the result must be correct because
##' none of the condition can be true by chance.
##'  
##' @title Test that a Provided Graident, Jacobian or Hessian is in
##'     Accordance with a Numeric Hessian
##'
##' @param der The provided vector, matrix or array of derivatives to
##'     be checked.
##' 
##' @param derNum The numeric vector, matrix or array of derivatives.
##' 
##' @param PREC A precision parameter. A value of \code{0.05} seems
##'     good in practice.
##' 
##' @return A logical which is \code{TRUE} if the two objects are
##'     close enough.
##'
##' @export
##' @keywords internal
testNumDeriv <- function(der, derNum, PREC = 0.05,
                          type = c("gradient", "hessian")) {
    if (any(is.na(der)) || any(is.na(derNum))) return(FALSE)
    type <- match.arg(type)
    small <- c("gradient" = 1e-3, "hessian" = 1e-2)[type] 
    (max(abs(der - derNum)) < small / PREC) ||
        (max(abs(der - derNum) / (abs(der) + 1e-9)) < small / PREC)
}

##' 
##' @title Check New Data for Prediction
##' 
##' @param object A \code{TSGEV} object.
##'
##' @param newdata An atomic vector or a data frame providing the data
##'     for a prediction.
##'
##' @param design Not used yet.
##' 
##' @return A data frame with the variables required for prediction.
##'     This must always contain a column with is name given by
##'     \code{object$date} and with class \code{"Date"}. If
##'     \code{object} contains \code{TSVars}, the corresponding
##'     variables must be provided as well, but with no check on their
##'     class.
##'
##' @keywords internal
##' 
## @examples
## df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
## fit <- TVGEV(data = df, response = "TXMax", date = "Date",
##              design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##              loc = ~ t1 + t1_1970)
## try(checkNewdata.TVGEV(fit, data.frame(date = "2026-01-20")))
## checkNewdata.TVGEV(fit, data.frame(Date = "2026-01-20"))
## 
## \dontrun{
##    data(fremantle)
##    df <- within(fremantle, Date <- as.Date(paste0(Year, "-01-01")))
##    fit <- TVGEV(data = df, response = "SeaLevel", date = "Date",
##                 design = polynomX(date = Date, degree = 1),
##                 loc = ~ t1 + SOI, trace = 2)
##    try(checkNewdata.TVGEV(fit, data.frame(Date = "2026-01-20", x = 2)))
##    checkNewdata.TVGEV(fit, data.frame(Date = "2026-01-20", SOI = 2.1))
## }
## 
checkNewdata.TVGEV <- function(object, newdata = NULL, design = FALSE) {

    
    if (length(object$TSVars)) {
        if (is.null(newdata)) {
            newdata <- object$data
            fDate <- object$fDate
        }
        newdata <- as.data.frame(newdata)
        requiredVars <- c(object$date, object$TSVars)
        notFound <- setdiff(requiredVars, names(newdata))
        if (length(notFound)) {
            stop("Required variables not found: ",
                 paste(notFound, collapse = ", "))
        }
    } else {
        if (is.null(newdata)) {
            newdata <- object$data[[object$date]]
            fDate <- object$fDate
        }
        if (is.atomic(newdata)) {
            newdata <- data.frame(as.Date(newdata))
            colnames(newdata) <- object$date
            fDate <- format(newdata[[object$date]])
        } else {
            newdata <- as.data.frame(newdata)
            if (ncol(newdata) == 1) {
                if (names(newdata) != object$date)
                    stop("Bad column name in 'newdata'.",
                         " Expecting '", object$date, "'.")
                newdata[[object$date]] <- as.Date(newdata[[object$date]])
            } else {
                if (is.na(m <- match(object$date, names(newdata))))
                    stop("Bad column names in 'newdata'.",
                         " Expecting to find '", object$date, "'.")
                newdata[[object$date]] <- as.Date(newdata[[object$date]])
                newdata <- newdata[, m, drop = FALSE]
            }
            fDate <- format(newdata[[object$date]])
        }
           
    }

    newdata
}
