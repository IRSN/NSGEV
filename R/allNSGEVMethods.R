##*****************************************************************************
##' Plot method for \code{NSGEV} objects.
##'
##' @title Plot method for \code{NSGEV} object
##'
##' @param x A NSGEV object
##'
##' @param y Not used.
##'
##' @param which An integer. The value \code{1} will lead to a plot
##' with the observation number used for the \eqn{x} axis. With the
##' value \code{2} the \eqn{x} axis is the covariate, assuming that
##' there is only one covariate.
##'
##' @param pLim A numeric vector with two probabilities used to set
##' the \code{ylim} parameter of \code{plot.default}. The model
##' quantiles for these probabilites are computed and added to the
##' range used to choose \code{ylim}.
##'
##' @param showQuant Not used yet.
##' 
##' @param ... Other arguments to be passed to methods.
##' 
##' @return Nothing.
##' 
##' @author Yves Deville
##'
##' @examples
##' example(NSGEV)
##' plot(ns1, which = 1)
##' plot(ns1, which = 2)
##' 
plot.NSGEV <- function(x, y, which = 1,
                       pLim = c(0.01, 0.99),
                       showQuant = c(0.5, 0.95, 0.99),
                       ...) {
    
    if (which == 1L) {

        n <- 60
        ysim <- simulate(x, nsim = n)
        
        if (!is.null(x$response)) {
            r <- range(x$response, ysim)
            ylab <- x$reponseName
            main <- "observed and simulated responses"
        }
        else {
            r <-range(x$response, ysim)
            ylab <- ""
            main <- sprintf("%d simulated responses", n)
        }
            
        matplot(ysim, type = "l", col = "gray", ylim = r,
                xlab = "obs.", ylab = ylab,
                main = main)
        
        if (!is.null(x$response)) {
            lines(x = x$response, type = "o", pch = 16)
        }
       
    } else if (which == 2L) {
        
        ## =============================================================
        ## Does it make sense with more than one predictor?
        ## XXX for the quantile function use a new dataset with a fine
        ## grid and sorted values
        ## =============================================================
        if (length(x$predNames) == 0L) {
            stop("'which = 2' is only possible when at least one ",
                 "predictor exists")
        }
        q0 <- quantile(x, probs = pLim)
        qLines <- (!is.logical(showQuant) || showQuant) && length(showQuant)
       
        if (!is.null(x$response)) {
            r <- range(x$response, q0)
            y <- x$response
            main <- "observed responses"
            ylab <- x$reponseName
        } else {
            r <-range(x$response, q0)
            y <- drop(simulate(x, nsim = 1))
            main <- "simulated responses"
            ylab <- x$reponseName
        }
        plot(x = x$data[ , x$predNames],  y = y, type = "p",
             ylim = r,
             pch = 21, col = "orangered", bg = "gold",
                xlab = x$predNames, main = main)
        if (qLines) {
            q <- quantile(x, probs = showQuant)
            matlines(x = x$data[ , x$predNames],  y = q, lwd = 2,
                     lty = 1:3, col = "darkgray")
        }
 
    }
}

##*****************************************************************************
##' Extract the vector of coefficients from a NSGEV object.
##'
##' @aliases coef.TVGEV
##' 
##' @title Coefficients of a \code{NSGEV} object
##'
##' @param object A \code{NSGEV} object.
##' 
##' @param type Character. For \code{type = "psi"}, the vector
##' \eqn{\boldsymbol{\psi}}{\psi} of model parameters is
##' returned. When instead \code{type} is \code{"theta"}, the matrix
##' of GEV parameters \eqn{\boldsymbol{\theta}_i}{\theta_i} is
##' returned, with one row by block (or observation) and one column
##' for each of the GEV parameters \code{"loc"}, \code{"scale"} and
##' \code{"shape"}.
##'
##' @param ... Not used yet.
##' 
##' @return Vector \eqn{\mathbf{\psi}}{\psi} of coefficients, or
##' matrix with the GEV parameters \eqn{\mathbf{\theta}_i}{\theta_i}
##' as its rows.
##' 
coef.NSGEV <- function(object, type = c("psi", "theta"),  ...) {
    type <- match.arg(type)
    if (type == "psi") return(object$estimate)
    else return(object$theta)
}

##*****************************************************************************
##' Names of the parameters of a model.
##'
##' @title Names of the Parameters of a Statistical Model
##'
##' @param object A parametric statistical model for which the
##' parameters must have names.
##'
##' @param ... Not used yet.
##'
##' @return The character vector of the names of the parameters of the
##' model.
##' 
parNames <- function(object, ...) {
   UseMethod("parNames")
}

##*****************************************************************************
##' Names of the Parameters of a Model.
##'
##' @title  Names of the parameters of a Statistical Model
##'
##' @param object A parametric statistical model for which the
##' parameters must have names.
##'
##' @param ... Not used yet.
##'
##' @return The character vector of the names of the parameters of the
##' model.
parNames.default <- function(object, ...) {

    if (is.list(object)) {
      nm <- names(object)
      ind <- (tolower(nm) == "parnames")
      if (sum(ind) == 1) return(object[[nm[ind]]])
   }
   
   invisible(NULL)
   
}

##*****************************************************************************
##' Summary method for \code{NSGEV} objects.
##'
##' @title Summary Method for \code{NSGEV} Objects
##'
##' @param object A \code{NSGEV} object.
##'
##' @param ... Not used yet.
##'
##' @return A list with the elements of \code{objects}
##' and some more that can be displayed when \code{summary}
##' is invoked.
summary.NSGEV <- function(object, ...) {
    
    res <- object
    formText <-
        paste(rep("    ", 3),
              paste(sprintf("%4s : ", names(object$formulas)),
                    object$formulas),
              collapse = "\n")
    res$formText <- formText
    class(res) <- "summary.NSGEV"
    res
}

##*****************************************************************************
print.summary.NSGEV <- function(x, ...) {
    
    cat("o Names of parameters (psi):\n   ",
        paste(sprintf("\"%s\"", x$parNames), collapse = ", "),
        "\n\n")
    if (length(x$predNames)) {
        cat("o Names of predictors (x):\n   ",
            paste(sprintf("\"%s\"", x$predNames), collapse = ", "),
            "\n\n")
    } else {
        cat("o No predictors (x)\n")
    }
    cat("o Formulas for GEV parameters:\n", x$formText, "\n\n")
    cat("o Coefficients:\n")
    print(x$estimate)
    cat("\n")

    cat("o GEV parameters:\n")
    d <- max(-ceiling(log(apply(x$theta, 2,
                           function(x) max(abs(x))) / 100, base = 10)))
    mat <- rbind(apply(x$theta, 2, range),
                 apply(x$theta, 2, mean))
    mat <- round(mat, digits = d)
    rownames(mat) <- c("min", "max", "mean")
    colnames(mat) <- sprintf("    %4s :", colnames(mat))
    print(t(mat))
    
}

##*****************************************************************************
print.NSGEV <- function(x, ...) {
    print(summary(x))
}
