
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

## ***************************************************************************
##' Bootstrap an object representing a fitted model.
##'
##' @title Bootstrap an Object
##' 
##' @param object An object representing a fitted model.
##'
##' @param ... Further arguments for methods.
##'
##' @return Bootstrap results
##' 
bs <- function(object, ...) {
    UseMethod("bs")
}

## ****************************************************************************
##' Marginal cumulative distribution functions.  
##'
##' @title Marginal Cumulative Distribution Functions
##'
##' @param x An object reprensenting a fitted model.
##'
##' @param ... Further arguments for methods.
##'
##' @return A structure containing discretized versions of the
##' marginal distribution functions for the model \code{x}.
##' 
cdf <- function(x, ...) {
    UseMethod("cdf")
}

## *************************************************************************
##' Compute the moments of the marginal distributions attached to a
##' model object.
##' 
##' @title Moments of Marginal Distributions
##' 
##' @param x An object representing a model.
##'
##' @param which The order of the moment.
##'
##' @param ... Further arguments for methods.
##'
##' @return A vector of moments.
##' 
moment <- function(x, which, ...) {
    UseMethod("moment")
}
