##' Compute vector(s) of GEV Parameters \code{theta} from vector(s)
##' of model parameters \code{psi}.
##' 
##' @title Compute GEV Parameters from Model Parameters
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
