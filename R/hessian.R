
## ****************************************************************************
##' Exact Hessian of the negative log-likelihood function of a
##' \code{TVGEV} object.
##'
##' The Hessian matrix computation relies on the exact second order
##' derivatives of the GEV density w.r.t the vector
##' \eqn{\boldsymbol{\theta}}{\theta} of GEV parameters, see
##' \code{\link[nieve]{dGEV}}. The derivatives w.r.t. the vector of model
##' parameters \eqn{\boldsymbol{\psi}}{\psi} are computed by chain rule
##' using the linear inverse link.
##' 
##' @title Hessian of the Negative Log-Likelihood of a \code{TVGEV} Object
##'
##' @param psi Value of the parameter vetor.
##'
##' @param object A \code{TVGEV} object.
##'
##' @return The Hessian matrix of the negative log-likelihood. When
##' \code{psi} is taken as the estimate of the model given in
##' \code{object}, the matrix should be positive semidefinite.
##' 
##' @note For now, only a fairly limited use is made of the Hessian,
##' hence we do not strive to optimise the computation. Note that the
##' Hessian used in the \code{\link{TVGEV}} estimation function
##' (creator of \code{TVGEV} objects) and providing the covariance
##' matrix \code{vcov} by inversion is a \emph{numeric} one for now.
##'
##' @section News: From version 0.1.8, the \code{negLogLikFun} function
##'     (closure) shipped with a \code{TVGEV} optionally comptes the
##'     Hessian.
##' 
##' @importFrom nieve dGEV
##'
##' @export
##' 
##' @examples
##' \dontrun{
##'    example(TVGEV)
##'    H1 <- d2negLogLikFun(res1$estimate, object = res1)
##'    H1 - res1$hessian
##'    H2 <- d2negLogLikFun(res2$estimate, object = res2)
##'    H2 - res2$hessian
##' }
d2negLogLikFun <- function(psi, object) {

    if (!inherits(object, "TVGEV")) {
        stop("'object' must inherit from the \"TVGEV\" class")
    }
    
    yVal <- object$data[object$indVal , object$response]
    
    parNames.GEV <- c("loc", "scale", "shape")
    nVal <- sum(object$indVal)
    thetaVal <- array(NA, dim = c(nVal, 3L),
                      dimnames = list(NULL, parNames.GEV))
    XVal <- list()
    
    for (nm in parNames.GEV) {
        if (!object$isCst[nm]) {
            XVal[[nm]] <- object$X[[nm]][object$indVal, , drop = FALSE]
        } else {
            XVal[[nm]] <- array(1.0, dim = c(nVal, 1L))
        }
        thetaVal[ , nm] <- XVal[[nm]] %*% psi[object$ind[[nm]]]
    }
    
    logL <-
        nieve::dGEV(yVal,
                    loc = thetaVal[ , "loc", drop = FALSE],
                    scale = thetaVal[ , "scale", drop = FALSE],
                    shape = thetaVal[ , "shape", drop = FALSE],
                    log = TRUE, deriv = TRUE, hessian = TRUE)
    
    H1 <- -attr(logL, "hessian")
    d2nL <- array(0.0, dim = c(object$p, object$p),
                  dimnames = list(object$parNames, object$parNames))
    
    for (j in 1L:3L) {

        nmj <- parNames.GEV[j]
        indj <- object$ind[[nmj]]
        
        for (k in j:3L) {

            nmk <- parNames.GEV[k]
            indk <- object$ind[[nmk]]
            s <- array(0.0, dim = c(object$pp[nmj], object$pp[nmk])) 
            
            for (i in 1:nVal) {
                s <- s - tcrossprod(XVal[[nmj]][i, ],
                                XVal[[nmk]][i, ]) *
                    attr(logL, "hessian")[i, nmj, nmk]
            }
            d2nL[indj, indk] <- d2nL[indk, indj] <- s
        }
        
        
    }
        
    return(d2nL)
    
}



