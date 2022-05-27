## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the implementation of the GEV distribution (C code used via
## .Call)
##
## The most difficult thing seems to be to have the good hessian of
## the quantile function when the shape 'xi' is very small, say 'xi' =
## 1e-4 or less. The results depends much on the value of 'eps' in the
## code 'GEV.c'. The value of eps should not be too small and a good
## choice seems to be around eps = 2e-4.
##
## NOT every thing is tested. We do not test that 'lower.tail' works
## as expected, ...
##
## ***************************************************************************

library(NSGEV)
library(numDeriv)
library(testthat)
context("GEV")

set.seed(1234)

## ==========================================================================
## check that the GEV quantile and distribution functions are consistent
## ==========================================================================

n <- 100
mu <- rnorm(1)
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)

x <- as.vector(rGEV(n = n, loc = mu, scale = sigma, shape = xi))
F <- pGEV(x, loc = mu, scale = sigma, shape = xi)
e <- x - qGEV(p = F, loc = mu, scale = sigma, shape = xi)

test_that(desc = "Consistency of the GEV cdf and quantile funs #1",
          expect_lt(max(abs(e)), 1e-10))

p <- runif(n)
q <- qGEV(p, loc = mu, scale = sigma, shape = xi)
e <- p - pGEV(q, loc = mu, scale = sigma, shape = xi)
test_that(desc = "Consistency of the GEV cdf and quantile funs #2",
          expect_lt(max(abs(e)), 1e-10))

## ==========================================================================
## check that the GEV density and distribution functions are
## consistent. Note that we do not check that the computation is 
## precise but rather that there is no error in the code.
## ==========================================================================

n <- 100
mu <- rnorm(1)
sigma <- rexp(1)

for (xi in c(-0.2, -0.1, -1e-4, -1e-7, 0.0, 1e-7, 1e-4, 0.1, 0.2)) {
    
    x <- rGEV(n, loc = mu, scale = sigma, shape = xi)
    
    F <- function(x) {
        pGEV(x, loc = mu, scale = sigma, shape = xi)
    }
    
    fval <- dGEV(x, loc = mu, scale = sigma, shape = xi)
    eps <- 1e-6
    e <- fval - (F(x + eps) - F(x - eps)) / 2 / eps
    
    test_that(desc = "Consistency of the GEV cdf density funs",
              expect_lt(max(abs(e)), 1e-6))
}


## ==========================================================================
## check the gradient and the Hessian of the log-density
## Note that the numerical error on the Hessian is larger than that on the
## gradient.
## ==========================================================================

funs <- list("log-density" = dGEV, "distribution" = pGEV, "quantile" = qGEV)
Hessian <- c("log-density" = TRUE, "distribution" = FALSE, "quantile" = TRUE)

n <- 1
mu <- rnorm(1)
sigma <- rexp(1)

for (nm in names(funs)) {

    Args <- ArgsVal <- as.list(formals(funs[[nm]]))
    
    if (nm == "log-density") {
        Args[["log"]] <- ArgsVal[["log"]] <- TRUE
    }
    if ("deriv" %in% names(Args)) ArgsVal[["deriv"]] <- TRUE
    if ("hessian" %in% names(Args)) {
        ArgsVal[["hessian"]] <- TRUE
        hessian <- TRUE
    } else hessian <- FALSE
    

    for (xi in c(-0.2, -0.1, -1e-4, -1e-7, 0.0, 1e-7, 1e-4, 0.1, 0.2)) {

        if (nm == "quantile") {
            x <- runif(1)
        } else {
            x <- as.vector(rGEV(n = 1, loc = mu, scale = sigma, shape = xi))
        }
        
        Args[[1]] <- ArgsVal[[1]] <- x
        
        theta0 <- c(loc = mu, scale = sigma, shape = xi)
        
        f <- function(theta) {
            Args[["loc"]] <- theta[1]
            Args[["scale"]] <- theta[2]
            Args[["shape"]] <- theta[3]
            do.call(funs[[nm]], Args)
        }

        ArgsVal[["loc"]] <- mu
        ArgsVal[["scale"]] <- sigma
        ArgsVal[["shape"]] <- xi

        fval <- do.call(funs[[nm]], ArgsVal)
        
        ## =====================================================================
        ## Check the gradient
        ## =====================================================================
        
        Jnum <- jacobian(func = f, x = theta0)
        J <- attr(fval, "gradient")
        
        cond <- (max(abs(J - Jnum)) < 4e-4) ||
            (max(abs(J - Jnum) / (abs(J) + 1e-9)) < 4e-3)
        
        test_that(desc = sprintf("Gradient of the GEV %s", nm),
                  expect_true(cond))

        ## =====================================================================
        ## Note that for the Hessian we do not expect a precision
        ## as high as for the gradient
        ## =====================================================================
        
        if (hessian) {
            Hnum <- hessian(func = f, x = theta0)
            H <- drop(attr(fval, "hessian"))
            cond <- (max(abs(H - Hnum)) < 4e-2) ||
                (max(abs(H - Hnum) / (abs(H) + 1e-9)) < 4e-2)
            if (is.na(cond) || !cond) {
                cat(sprintf("testing hessian xi = %6.3f\n", xi))
                print(H)
                print(Hnum)
            }
            test_that(desc = sprintf("Hessian of the GEV %s", nm),
                      expect_true(cond))
        }
            
    }
}
    
