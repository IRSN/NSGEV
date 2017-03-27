
reshapeGEV <- function (x, loc, scale, shape, matrix = FALSE) {
   
   n <- length(x)
   
   if (n > 1L) {
      msg <- paste("when 'x' has length n > 1, 'loc', 'scale' and 'shape'",
                   " must have length 'n' or 1")
      if (length(loc) != n) {
         if (length(loc) != 1L) stop(msg)
         loc <- rep(loc, n)
      }
      if (length(scale) != n) {
         if (length(scale) != 1L) stop(msg)
         scale <- rep(scale, n)
      }
      if (length(shape) != n) {
         if (length(shape) != 1L) stop(msg)
         shape <- rep(shape, n)
      }
      nms <- names(x)
   } else {

       n <- max(c(length(loc), length(scale), length(shape)))
       
       if (n > 1L) {
          
          if (((length(loc) > 1L) && (length(loc) < n)) ||
              ((length(scale) > 1L) && (length(scale) < n)) ||
              ((length(shape) > 1L) && (length(shape) < n)) ) {
              stop("when 'x' has length 1, the lengths of the non-scalar elements",
                   "among 'loc', 'scale' and 'shape' must be the same")
          }
          x <- rep(x, n)
          loc <- rep(loc, length.out = n)
          scale <- rep(scale, length.out = n)
          shape <- rep(shape, length.out = n) 
      }
       
   }
   
   if (matrix)  {

       return(cbind("x" = x, "loc" = loc, "scale" = scale, "shape" = shape))
      
   } else {

       if (length(x) > 1L) nms <- names(x)
       else if (length(loc) > 1L) nms <- names(loc)
       else if (length(scale) > 1L) nms <- names(scale)
       else if (length(shape) > 1L) nms <- names(shape)
       else nms <- ""
       
       return(list("n" = n, "x" = x,
                   "loc" = loc, "scale" = scale, "shape" = shape,
                   "nms" = nms))
   }
}

##' Density, distribution function, quantile function and random
##' generation for the Generalized Extreme Value (GEV) distribution
##' with parameters \code{loc}, \code{scale} and \code{shape}.
##'
##' @name GEV
##' @rdname GEV
##' 
##' @title Density, Distribution Function, Quantile Function and
##' Random Generation for the Generalized Extreme Value (GEV)
##' Distribution
##'
##' 
##' @param loc Location parameter. Numeric vector with suitable
##' length, see \bold{Details}.
##'
##' @param scale Scale parameter. Numeric vector with suitable length,
##' see \bold{Details}.
##'
##' @param shape Shape parameter. Numeric vector with suitable length,
##' see \bold{Details}.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##' returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##' computed value w.r.t. the parameter vector is computed.
##'
##' @param impl Character. Choose between a R implementation and a C
##' implementation. The second one should be faster but less checks
##' have been done on it for now.
##'
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
##' P[X <= x], otherwise, P[X > x].
##'
##' @return A numeric vector with length equal to the length of the
##' first argument or of the parameters. When \code{deriv} is
##' \code{TRUE}, the returned value has an attribute named
##' \code{"gradient"} which is a matrix with \eqn{n} lines and \eqn{3}
##' columns containing the derivatives. A row contains the partial
##' derivatives of the corresponding element w.r.t. the three
##' parameters \code{loc} \code{scale} and \code{shape} in that order.
##'
##' @details For the \code{d}, \code{p} and \code{q} functions, let
##' \code{n} be the length of the first element. If \code{n > 1}, then
##' each of the parameters \code{loc}, \code{scale} and \code{shape}
##' must be of length \code{1} or \code{n}: in the first case, it will
##' be recycled to have length \code{n}. When \code{n = 1} the largest
##' length met for the three parameters, say \code{np} is used. The
##' two other parameters must then be of length \code{1} or \code{np},
##' and are given the length \code{np}, as well as the first argument.
##' Note that only vectors of length one are actually recycled.
##'
##' @author Yves Deville
##'
##' @examples
##' ti <- 1:10; names(ti) <- 2000 + ti
##' mu <- 1.0 + 0.1 * ti
##' ## simulate 40 paths
##' y <- rGEV(n = 40, loc = mu, scale = 1, shape = 0.05)
##' matplot(ti, y, type = "l", col = "gray")
##' lines(ti, apply(y, 1, mean))
dGEV <- function(x, loc = 0.0, scale = 1.0, shape = 0.0, log = FALSE,
                 deriv = FALSE, impl = c("C", "R")) {

    impl <- match.arg(impl)
    
    if (impl == "C") {
        
        res <- .Call("Call_dGEV",
                     as.double(x),
                     as.double(loc),
                     as.double(scale),
                     as.double(shape),
                     as.integer(log),
                     as.integer(deriv),
                     PACKAGE = "NSGEV")

        n <- length(res)
        if (deriv) {
            attr(res, "gradient") <-
                array(attr(res, "gradient"),
                      dim = c(n, 3),
                      dimnames = list(rownames(x), c("loc", "scale", "shape")))
        }
        return(res)
        
    } 
    
    L <- reshapeGEV(x = x, loc = loc, scale = scale, shape = shape,
                    matrix = TRUE)

    n <- nrow(L)
    d <- rep(NA, n)
    z <- (L[ , "x"] - L[ , "loc"]) / L[ , "scale"]
    
    if (deriv) {
        grad <- array(NA, dim = c(n, 3L),
                      dimnames = list(rownames(L), c("loc", "scale", "shape")))
    }

    nax <- is.na(L[ , "x"])
    
    ## Gumbel xi = 0.0
    ind <- (!nax  & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) < 1e-6))
    
    if (any(ind)) {
        z_ind <- z[ind]
        scale_ind <- L[ind, "scale"]
        emz_ind <- exp(-z_ind)
        d[ind] <- -log(scale_ind) - z_ind - emz_ind
        if (deriv) {

            ## grad[ind, ] <- c("loc" = (1.0 - emz_ind) / scale_ind,
            ##                  "scale" = (-1.0 + z_ind * (1.0 - emz_ind)) / scale_ind,
            ##                  "shape" = z_ind * z_ind* (1 - emz_ind) / 2.0 - z_ind)
            grad[ind, "loc"] <- (1.0 - emz_ind) / scale_ind
            grad[ind, "scale"] <- (-1.0 + z_ind * (1.0 - emz_ind)) / scale_ind
            grad[ind, "shape"] <- z_ind * z_ind* (1 - emz_ind) / 2.0 - z_ind
        
        }
    }
    ## non-Gumbel xi != 0.0
    ind <- (!nax  & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) >= 1e-6))
    
    if (any(ind)) {
        d_ind <- rep(-Inf, sum(ind))
        if (deriv) {
            grad_ind <- array(0, dim = c(sum(ind), 3L),
                              dimnames = list(rownames(L), c("loc", "scale", "shape")))
        }
        z_ind <- z[ind]
        V_ind <- 1.0 + L[ind, "shape"] * z_ind
        xi_ind <- L[ind, "shape"]
        sigma_ind <- L[ind, "scale"]
        ind2 <- (V_ind > 0)
        if (any(ind2)) {
            d_ind[ind2] <- -log(sigma_ind[ind2]) -
                V_ind[ind2]^(-1.0 / xi_ind[ind2]) - 
                    (1.0 / xi_ind[ind2] + 1.0) * log(V_ind[ind2])
            if (deriv) {
                W_ind <- V_ind^(-1.0 / xi_ind)
                U_ind <- (1.0 + xi_ind - W_ind) / V_ind / sigma_ind
                ## grad_ind[ind2, ] <-
                ##     c("loc" =  U_ind[ind2],
                ##       "scale" = -1.0 / sigma_ind[ind2] + z_ind[ind2] * U_ind[ind2],
                ##       "shape" =  log(V_ind[ind2]) * (1.0 - W_ind[ind2]) /
                ##           xi_ind[ind2] / xi_ind[ind2] -
                ##           z_ind[ind2] * U_ind[ind2] * sigma_ind[ind2] / xi_ind[ind2])
                grad_ind[ind2, "loc"] <- U_ind[ind2]
                grad_ind[ind2, "scale"] <- -1.0 / sigma_ind[ind2] + z_ind[ind2] * U_ind[ind2]
                grad_ind[ind2, "shape"] <- log(V_ind[ind2]) * (1.0 - W_ind[ind2]) /
                    xi_ind[ind2] / xi_ind[ind2] -
                        z_ind[ind2] * U_ind[ind2] * sigma_ind[ind2] / xi_ind[ind2]

                
            }
        }
        d[ind] <- d_ind
        if (deriv) {
            grad[ind] <- grad_ind
        }
    }
    
    if (!log) {
        d <- exp(d)
        if (deriv) {
            ## multiply the gradients by the density
            grad <- sweep(x = grad, MARGIN = 1L, STATS = d, FUN = "*")
        }
    }
    if (deriv) {
        attr(d, "gradient") <- grad
    }
    d
    
}

##' @rdname GEV
pGEV <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  deriv = FALSE) {


    L <- reshapeGEV(x = q, loc = loc, scale = scale, shape = shape,
                    matrix = TRUE)
    n <- nrow(L)
    
    if (deriv) {
        grad <- array(NA, dim = c(n, 3L),
                      dimnames = list(rownames(L), c("loc", "scale", "shape")))
    }
    
    p <- rep(NA, n)
    z <- (L[ , "x"] - L[ , "loc"]) / L[ , "scale"]
    nax <- is.na(L[ , "x"])
                 
    ## Gumbel xi = 0.0
    ind <- (!nax & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) < 1e-6))
    if (any(ind)) {
        if (deriv) {
            z_ind <- z[ind]
            emz_ind <- exp(-z_ind)
            p[ind] <-  exp(-emz_ind)
            Z_ind <-  emz_ind * p[ind]
            grad[ind, ] <- c("loc" = -Z_ind / L[ind, "scale"],
                             "scale" = -z_ind * Z_ind / L[ind, "scale"],
                             "shape" = -z_ind * z_ind * Z_ind / 2)
        } else {
             p[ind] <-  exp(-exp(-z[ind]))
        }
    }
    ## non-Gumbel xi != 0.0
    ind <- (!nax & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) >= 1e-6))
    if (any(ind)) {
        if (deriv) {
            grad_ind <- array(0, dim = c(sum(ind), 3L)) 
        }
        xi_ind <- L[ind, "shape"]
        ## set the value for V_ind <= 0.0. When xi > 0 this is 0.0, and for
        ## xi < 0.0 this is 1.0
        p_ind <- (xi_ind < 0.0)   
        V_ind <- 1.0 + L[ind, "shape"] * z[ind]
        ind2 <- (V_ind > 0.0)
        if (any(ind2)) {
            if (deriv) {
                xi_ind2 <- xi_ind[ind2]
                sigma_ind2 <- L[ind, "scale"][ind2]
                z_ind2 <- z[ind][ind2]
                V_ind2 <- V_ind[ind2]
                W_ind2 <- V_ind[ind2]^(-1.0 / xi_ind2)
                Z_ind2 <- W_ind2 * exp(-W_ind2)
                p_ind[ind2] <- exp(-W_ind2)
                grad_ind[ind2, ] <-
                    c("loc" =  -Z_ind2 / V_ind2 / sigma_ind2,
                      "scale" = -z_ind2 * Z_ind2 / V_ind2 / sigma_ind2,
                      "shape" = -Z_ind2 * (log(V_ind2) / xi_ind2  - z_ind2 / V_ind2) / xi_ind2)      
            } else {
                p_ind[ind2] <- exp(-V_ind[ind2]^(-1.0 / xi_ind[ind2]))
            }
        }
        p[ind] <- p_ind
        if (deriv) {
            grad[ind, ] <- grad_ind 
        }
    }
    if (!lower.tail)  {
        p <- 1 - p
        if (deriv) {
            grad <- -grad
        } 
    }
    if (deriv) {
        attr(p, "gradient") <- grad
    }
    p
}

##' @rdname GEV
qGEV <- function (p, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE,
                  deriv = FALSE) {
    
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
    
    L <- reshapeGEV(x = p, loc = loc, scale = scale, shape = shape,
                    matrix = TRUE)
    n <- nrow(L)
    
    if (deriv) {
        grad <- array(NA, dim = c(n, 3L),
                      dimnames = list(rownames(L), c("loc", "scale", "shape")))
    }
    q <- rep(NA, n)
    if (!lower.tail) L[ , "x"] <- 1.0 - L[ , "x"]

    nap <- is.na(p)
    
    ## Gumbel xi = 0.0
    ind <- (!nap & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) < 1e-6))
    if (any(ind)) {
        A_ind <- -log(L[ind, "x"])
        q[ind] <- L[ind, "loc"] - L[ind, "scale"] * log(A_ind)
        if (deriv) {
            logA_ind <- log(A_ind)
            grad[ind, ] <- c("loc" = rep(1.0, sum(ind)),
                             "scale" = -logA_ind,
                             "shape" = L[ind, "scale"] * logA_ind^2 / 2.0)
        }
    }
    ## non-Gumbel xi != 0.0
    ind <- (!nap & (L[ , "scale"] > 0.0) & (abs(L[ , "shape"]) >= 1e-6))
    if (any(ind)) {
        A_ind <- -log(L[ind, "x"])
        xi_ind <- L[ind, "shape"]
        V_ind <- (1.0 - A_ind^(-xi_ind)) / xi_ind
        q[ind] <- L[ind, "loc"] + L[ind, "scale"] * (A_ind^(-xi_ind) - 1.0) / xi_ind
        if (deriv) {
            grad[ind, ] <-
                c("loc" = rep(1.0, sum(ind)),
                  "scale" = -V_ind,
                  "shape" = L[ind, "scale"] * (V_ind - log(A_ind) * (-xi_ind * V_ind + 1.0)) / xi_ind)
        }
    }
    if (deriv) {
        attr(q, "gradient") <- grad
    }
    q
   
}

##' @rdname GEV
rGEV <- function (n, loc = 0.0, scale = 1.0, shape = 0.0) {
    
    if (any(is.na(loc)) || !all(is.finite(loc))) {
        stop("'loc' must contain non-NA finite numeric values")  
    }
    if (any(is.na(scale)) || any(scale <= 0) || !all(is.finite(scale))) {
        stop("'scale' must contain non-NA finite and positive numeric values")  
    }
    if (any(is.na(shape)) || !all(is.finite(shape))) {
        stop("'shape' must contain non-NA finite numeric values")  
    }
    
    L <- reshapeGEV(x = 1.0, loc = loc, scale = scale, shape = shape)
    
    r <- array(NA, dim = c(L$n, n),
               dimnames = list(L$nms, paste("sim", 1:n, sep = "")) )
    
    loc <- array(L$loc, dim = c(L$n, n))
    scale <- array(L$scale, dim = c(L$n, n))
    shape <- array(L$shape, dim = c(L$n, n))
    
    ## Gumbel xi = 0.0
    ind <- (abs(L$shape) < 1e-6)
    if (any(ind)) {
        nl <- sum(ind)
        n_ind <- nl * n
        U <- array(runif(n_ind), dim = c(nl, n))
        r[ind, ] <- loc[ind, drop = FALSE] - scale[ind, drop = FALSE] * log(U)
    }
    
    ## non-Gumbel xi != 0.0
    ind <- (abs(L$shape) >= 1e-6)
    if (any(ind)) {
        nl <- sum(ind)
        n_ind <- nl * n
        U <- array(runif(n_ind), dim = c(nl, n))
        U <- -log(U)
        r[ind, ] <- loc[ind, drop = FALSE] +
            scale[ind, drop = FALSE] * (U^(-shape[ind, , drop = FALSE]) - 1.0) /
                shape[ind, drop = FALSE]
    }
    r
}

