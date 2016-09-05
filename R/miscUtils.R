##' Translate a standard formula into a NSGEV formula.
##'
##' @title Translate a formula into a NSGEV formula
##' 
##' @param formula The formula to be tranlated.
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
