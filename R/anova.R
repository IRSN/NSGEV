##' Compute an analysis of deviance table for two nested \code{TVGEV}
##' objects.
##'
##' 
##' @title Analysis of Deviance Table for two Nested \code{TVGEV}
##'     Objects
##' 
##' @param object A \code{TVGEV} object as fitted by
##'     \code{\link{TVGEV}}.
##' 
##' @param object1 A \code{TVGEV} object such that \code{object} is
##'     nested in \code{object1}.
##' 
##' @param trace Level of verbosity. The value \code{0} prints
##'     nothing.
##' 
##' @param ... Not used yet.
##'
##' @return An object of class \code{"anova"} inheriting from class
##'     \code{"data.frame"}.
##'
##' @note The deviance of the models can not be interpreted: only the
##'     difference of the deviance is used.
##'
##' @section Caution: The distribution of the test statistic
##'     (difference between two deviances) is obtained on the basis of
##'     the \emph{large sample theory} which may not be applicable.
##'
##' @importFrom stats anova pchisq
##' @method anova TVGEV
##' @export
##' 
##' @examples
##'
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##'
##' ## fit without a break in trend
##' res0 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'                design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'                loc = ~ t1,
##'                estim = "nloptr")
##' 
##' ## the same but with a break
##' res1 <- TVGEV(data = df, response = "TXMax", date = "Date",
##'               design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'               loc = ~ t1 + t1_1970,
##'               estim = "nloptr")
##' anova(res0, res1)
##' 
anova.TVGEV <- function(object, object1, trace = 1L, ...) {
    
    if (missing(object) || missing(object1))
        stop("two models must be specified in 'object' and 'object1'")
    
    model0 <- deparse(substitute(object))
    model1 <- deparse(substitute(object1))
    models <- list(model0, model1)
    
    if (trace) {
 
    }
    
    ## check that the models have the same data. Could be improved
    ## by adding more tests (data values, ...)
    for (elt in list("nobs")) {
        if (object[[elt]] != object1[[elt]]) 
            stop("'object' and 'object1' must have the same '",
                 paste(elt, collapse = "$"), "' element") 
    }
     
    
    df <- c(object$df, object1$df)
    dev <- -2 * c(object$logLik, object1$logLik)
    dfDiff <- diff(df)
    devDiff <- diff(dev)
                  
    ## checks df
    if (dfDiff <= 0) stop("non-nested models. Bad order for models?")
         
    w <- -devDiff

    pVal <- pchisq(w, df = 1, lower.tail = FALSE) 
    
    table <- data.frame(df, dev,  c(NA, w), c(NA, pVal))
    dimnames(table) <- list(models, c("df", "deviance",  "W", "Pr(>W)"))
    
    structure(table, heading = c("Analysis of Deviance Table\n"),
              class = c("anova", "data.frame"))
    
    ## list(W = w, p.value = pVal, method = method)
    
}


