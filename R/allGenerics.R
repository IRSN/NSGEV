
## *************************************************************************
##' Compute vector(s) of GEV Parameters \code{theta} from vector(s)
##' of model parameters \code{psi}.
##' 
##' @title Compute the Matrix of GEV Parameters from the Vector of
##' Model Parameters 'psi' and the Data
##' 
##' @param model An object representing a model with GEV margins.
##'
##' @param psi A vector of parameters for the model.
##'
##' @param ... Not used yet.
##'
##' @return A vector or matrix of GEV parameters.
psi2theta <- function(model, psi, ...) {
    UseMethod("psi2theta")
}


bs <- function(object, ...) {
    UseMethod("bs")
}
