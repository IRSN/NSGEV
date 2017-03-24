##' Translate a standard formula into a NSGEV formula.
##'
##' @title Translate a formula into a NSGEV formula
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
##' @examples
##' 
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

##' Generates random names.
##'
##' This function is a utilisty function for tests concerning the
##' parsing of formulas.
##' 
##' @title Generates random names
##'
##' @param n Number of names.
##'
##' @param nchar Maximal number of characters for each name.
##' 
##' @return A character vector of valid names.
##'
##' @examples
##' set.seed(31415)
##' rnames(4)
##' rnames(4, nchar = 6)
rnames <- function(n, nchar = 3L) {
    df <- expand.grid(LETTERS, LETTERS, LETTERS)
    LET <- apply(X = df[sample(1:676, n) ,],
                 MARGIN = 1L,
                 FUN = function(x) x[sample(1L:3L, sample(1L:3L))])
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
