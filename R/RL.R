
## *************************************************************************
##' Return Levels for \code{NSGEV} and \code{TVGEV} Objects.
##'
##' In all cases, the GEV parameters
##' \eqn{\boldsymbol{\theta}_b}{\theta_b} corresponding to the model
##' and the rows \eqn{b} of \code{data} are computed. When
##' \code{sampleData} is \code{TRUE} the rows are assumed to form a
##' sample of the target distribution of the GEV parameters and will
##' be re-sampled, so the result will change in successive calls. When
##' instead \code{sampleData} is \code{FALSE}, the rows of \code{data}
##' are considered as the values of the covariates that will be
##' observed in the next 'future' blocks. Then the Return period can
##' not be greater than the number of rows in \code{'data'}.
##'
##' \itemize{
##'
##'    \item{When \code{type} is \code{"average"}, the RL is computed
##' for each value of the GEV parameter and the result is simply the
##' mean of the computed RLs.}
##' 
##'   \item{When \code{type} is \code{"exceed"} the RL is the value
##' \eqn{\rho} for which the number of exeedances over \eqn{\rho} has
##' unit expectation, as proposed by Parey et al. This value is found
##' by computing the expectation and solving a non-linear equation in
##' \eqn{\rho} with \code{uniroot} function.}
##' 
##' }
##' 
##' @title Return Levels for \code{NSGEV} and \code{TVGEV} objects
##'
##' @param period The Return period. Should be understood as an integer
##' multiple of the block duration \code{w}, and will be rounded else.
##' In some cases \code{period} must not exceed the number of rows in
##' \code{data}, see \bold{Details}.
##' 
##' @param model A \code{NSGEV} or \code{TVGEV} object.
##'
##' @param data To be used only when \code{model} has class
##' \code{"NSGEV"}. A data frame containing covariates for the model.
##' By default, the data frame attached to \code{model} is used.
##'
##' @param date To be used only when \code{model} has class
##' \code{"TVGEV"}. A vector that can be coerced to the \code{"Date"}
##' class.  By default, the date attached to \code{model} is
##' used.
##' 
##' @param psi Vector of coefficients. By default, the vector of
##' coefficients of the model is used.
##'
##' @param w Numeric scalar: duration of the blocks.
##'
##' @param type The type of RL wanted, see \bold{Details}.
##'
##' @param sampleData Logical. Only used for \code{NSGEV} objects. If
##' \code{TRUE}, the rows of \code{data} are assumed to form a sample
##' of the distribution of the covariates and will be re-sampled.
##'
##' @param deriv Logical. If \code{TRUE} the gradient of the RL w.r.t
##' the model parameters is returned.
##' 
##' @param rhoMin Minimal value for the evaluation of the function in
##' the zero-finding step when \code{type} is \code{"exceed"}. This
##' should be used if the zero is not found with the default value,
##' which can be diagnosed then with \code{plot = TRUE}.
##' 
##' @param rhoMax Maximal value for the evaluation of the function in
##' the zero-finding step when \code{type} is \code{"exceed"}. This
##' should be used if the zero is not found with the default value,
##' which can be diagnosed then with \code{plot = TRUE}.
##' 
##' @param plot Logical. If \code{type} is \code{"expect"}, a plot is
##' shown to check the results of the zero-finding step.
##'
##' @return The numeric value of the Return Level.
##'
##' @author Yves Deville
##'
##' %% note Consider the special case of a model simply defining a
##' %% linear trend time for the GEV location parameter \eqn{\mu} with a
##' %% constant scale \eqn{\sigma} and a constant shape \eqn{\xi}. If the
##' %% observations in \code{data} are the expectations for future blocks
##' %% conditional on past blocks, the computed Return Level.
##'
##' @section Caution: For now, \code{period} can only be a numeric
##' vector with length 1.
##'
##' @references
##'
##' Parey S., Hoang T.T.H., Dacunha-Castelle D. (2007) "Different ways
##' to compute temperature return levels in the climate change context".
##' \emph{Environmetrics}, \bold{21}, pp. 698-718.
##' 
##' @examples
##' ## =================
##' ## NSGEV examples
##' ## =================
##' example(NSGEV)
##' RL(ns1, period = 10, type = "average")
##' RL(ns1, period = 10, type = "exceed")
##' ## with derivative
##' RL(ns1, period = 10, type = "average", deriv = TRUE)
##' RL(ns1, period = 10, type = "exceed", deriv = TRUE)
##' 
##' ## check the zero-finding step
##' RL(ns1, period = 10, type = "exceed", deriv = TRUE, plot = TRUE)
##'
##' ## =================
##' ## TVGEV examples
##' ## ================
##' example(TVGEV)
##' RL(res1, period = 30)
##' nd <-  seq(from = as.Date("2000-01-01"), to = as.Date("2300-01-01"),
##'             by = "year")
##' RLe <- RL(res1, period = 200, date = nd, plot = TRUE)
##' RLa <- RL(res1, period = 200, date = nd, plot = TRUE, type = "average")
##' ## check the value of 'RLA'
##' q <- quantile(res2, prob = 1.0 - 1.0 / 200, date = nd)
##' plot(q)
##' mean(q[1:200])
##' 
RL <- function(model,
               period,
               data = NULL,
               date = NULL,
               psi = NULL,
               w = 1.0,
               type = c("exceed", "average"),
               sampleData = FALSE,
               deriv = FALSE,
               rhoMin = NULL,
               rhoMax = NULL,
               plot = FALSE) {

    type <- match.arg(type)
    p <- model$p
    
    if (period <= 1) stop("'period' must be > 1")
    
    if (is.null(psi)) psi <- model$estimate
    
    if (is(model, "NSGEV")) {
        
        if (!is.null(date)) {
            stop("when 'model' has class \"NSGEV\", the 'date' argument must\n",
                 "not be given. Use 'data' with a suitable data.frame instead.")
        }
        
        theta <- psi2theta(model = model, psi = psi, data = data, deriv = deriv,
                           checkNames = FALSE)
        
    } else if (is(model, "TVGEV")) {

        if (!is.null(data)) {
            stop("when 'model' has class \"TVGEV\", the 'data' argument must\n",
                 "not be given. Use 'date' with a suitable vector instead.")
        }
        if (sampleData) {
            stop("when 'model' has class \"TVGEV\", the 'sampleData' argument must\n",
                 "not be given.")
        }
        if (is.null(date)) date <- model$fDate
        date <- as.Date(date)
        theta <- psi2theta(model = model, psi = psi, date = date, deriv = deriv,
                           checkNames = FALSE)
    }

    nd <- nrow(theta)
    
    ## SOG: Save Our Gradient
    if (deriv) dtheta_dpsi <- attr(theta, "gradient")
    
    B <- period %/% w
    
    ## ========================================================================
    ## Rather than resampling the rows of theta, we use vectors of
    ## indices and weights
    ## =========================================================================

    if (!sampleData) {
        
        if (B > nd) {
            stop("When 'sampleData' is FALSE, the return period in 'period' ",
                 "should not be greater  than the number of rows in 'data'")
        } else {
            ## is this useful ????
            ## warning("Since 'sampleData' is FALSE, only the ", B, " first ",
            ##         " rows of 'data' are used")
            i_ind <- 1L:B
            n_ind <- B
            w_ind <- rep(1.0, B)
        }
        
    } else {
        ind <- sample(x = 1L:nd, size = B, replace = TRUE)
        t_ind <- table(ind)
        i_ind <- as.integer(names(t_ind))
        n_ind <- length(i_ind)
        w_ind <- as.numeric(t_ind) 
    }
    
    theta <- theta[i_ind, , drop = FALSE]
    if (deriv) dtheta_dpsi <- dtheta_dpsi[i_ind, , , drop = FALSE] 
    
    if (type == "average") {
        
        ##=====================================================================
        ## When type is "average", the RL 'rho' is simply the mean of
        ## the RLs corresponding to the rows of 'data', with possible
        ## resampling.
        ## =====================================================================

        q <- qGEV(1 - 1 / B, loc = theta[ , 1L], scale = theta[ , 2L],
                  shape = theta[ , 3L], deriv = deriv)
        
        rho <- weighted.mean(x = q, w = w_ind)

        if (deriv) {
            ## chain rule: dq / dpsi = (dq / dtheta) %*% (dtheta / dpsi)
            drho_dpsi <- array(0, dim = c(1, p))
            ## XXX to be improved by an array operation
            for (i in 1L:n_ind) {
                drho_dpsi <- drho_dpsi + 
                    tcrossprod(attr(q, "gradient")[i, ], 
                               dtheta_dpsi[i , , ]) * w_ind[i] / B
            }
            attr(rho, "gradient") <- drho_dpsi 
        }

        return(rho)
        
    } else if (type == "exceed") {
        
        ##=====================================================================
        ##  When type is "exceed", the RL 'rho' is found by solving a
        ##  non-linear equation.
        ## =====================================================================
        
        g <- function(rho) {
            F <- pGEV(rho, loc = theta[ , 1L], scale = theta[ , 2L],
                      shape = theta[ , 3L])
            s <- sum(F * w_ind) - (B - 1.0)
        }
        
        if (is.null(rhoMin) || is.null(rhoMax)) {
            ## to find the initial interval
            q <- qGEV(1 - 1 / B, loc = theta[ , 1L], scale = theta[ , 2L],
                      shape = theta[ , 3L], deriv = FALSE)
            if (is.null(rhoMin)) rhoMin <- min(q)
            if (is.null(rhoMax)) rhoMax <- max(q)

            ## added on 2017-08-18. This tyîcally occurs when the
            ## parameter falls outside of the admissible region, e.g.
            ## when the scale is negative.
            if (is.na(rhoMin) || is.na(rhoMax)) {
                rho <- NA
                if (deriv) {
                    attr(rho, "gradient") <- rep(NA, model$p)
                }
                return(rho)
            }
        }

        ## cat("XXX\n")
        ## print(c(g(rhoMin), g(rhoMax)))
        
        res <- try(uniroot(f = g, interval = c(rhoMin, rhoMax)))
        
        if (plot) {
            rhoCand <-
                exp(seq(from = log(1.1), to = log(rhoMax), length.out = 100))
            valCand <- sapply(rhoCand, g)
            
            plot(rhoCand, valCand, type = "o", pch = 16, cex = 0.5,
                 lwd = 2, col = "SteelBlue3",
                 main = "zero-finding", xlab = "rho", ylab = "g")
            abline(h = 0, col = "orangered")
            if (class(res) != "try-error") {
                abline(v = res$root, col = "orangered")
            }
            mtext(side = 1, at = res$root,
                  line = 1.2, text = "rho", col = "orangered")
        }
        
        if (class(res) == "try-error") {
            print(c(rhoMin, rhoMax))
            stop("XXXX")
        }
             
        if (abs(res$f.root) > 1e-3) stop("no solution found")
        rho <- res$root
            
        ## =====================================================================
        ## It the gradient is needed we use implicit function
        ## derivation and chain rule
        ## =====================================================================
        if (deriv) {
            
            F <- pGEV(rho, loc = theta[ , 1L], scale = theta[ , 2L],
                      shape = theta[ , 3L], deriv = TRUE)

            ## dg_rho is a scalar (obtained by summing on blocks)
            dg_drho <- sum(dGEV(rho, loc = theta[ , 1L], scale = theta[ , 2L],
                                shape = theta[ , 3L], deriv = FALSE) * w_ind)
            if (plot) {
                ## show the derivative
                abline(a = -rho * dg_drho, b = dg_drho, col = "SpringGreen4")
            }
            
            ## chain rule: dF / dpsi = (dF / dtheta) %*% (dtheta / dpsi)
            dg_dpsi <- array(0, dim = c(1, p))
            ## XXX to be improved by an array operation
            for (i in 1L:n_ind) {
                dg_dpsi <- dg_dpsi + 
                    tcrossprod(attr(F, "gradient")[i, ], 
                               dtheta_dpsi[i , , ]) * w_ind[i] 
            }
            attr(rho, "gradient") <- -dg_dpsi / dg_drho 
        }
        
        rho
        
    }
    
}

