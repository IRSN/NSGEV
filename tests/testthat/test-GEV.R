context("GEV")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the implementation of the GEV distribution (C code used via
## .Call)
## ***************************************************************************

library(numDeriv)
library(testthat)

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

for (xi in c(-0.1, -1e-4, 0.0, 1e-4, 0.1)) {
    
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

n <- 1
mu <- rnorm(1)
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)
x <- as.vector(rGEV(n = n, loc = mu, scale = sigma, shape = xi))

f <- function(theta) {
    dGEV(x, loc = theta[1], scale = theta[2], shape = theta[3], log = TRUE)
}

fval <- dGEV(x = x, loc = mu, scale = sigma, shape = xi, log = TRUE,
             deriv = TRUE, hessian = TRUE)

theta0 <- c(loc = mu, scale = sigma, shape = xi)
Jnum <- jacobian(func =f, x = theta0)
Hnum <- hessian(func = f, x = theta0)

test_that(desc = "Gradient of the GEV log-density",
            expect_lt(max(abs(attr(fval, "gradient") - Jnum)), 1e-6))

test_that(desc = "Hessian of the GEV log-lik",
            expect_lt(max(abs(drop(attr(fval, "hessian")) - Hnum)),
                      1e-4))

## ==========================================================================
## check that the Hessian of the log-likelihood is OK
## ==========================================================================

n <- 300
mu <- rnorm(1)
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)
y <- as.vector(rGEV(n = n, loc = mu, scale = sigma, shape = xi))
theta <- c(loc = mu, scale = sigma, shape = xi)

## this is for numerical differentiation
f2 <- function(theta, y) {
    sum(dGEV(x = y, loc = theta[1], scale = theta[2],
             shape = theta[3], log = TRUE))
}

## this is to get the Hessian as an attribute
fval <- dGEV(x = y, loc = mu, scale = sigma, shape = xi, log = TRUE,
             deriv = TRUE, hessian = TRUE)

H2 <- apply(attr(fval, "hessian"), MARGIN = c(2, 3), FUN = sum)
H2num <- hessian(func = f2, x = theta, y = y)

test_that(desc = "Hessian of the GEV log-lik",
            expect_lt(max(abs(H2 - H2num)), 1e-3))


## ==========================================================================
## Check the gradient and the Hessian of the quantile function
## ==========================================================================

n <- 1
mu <- rnorm(1,)
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)
p <- runif(n)

f <- function(theta) {
    qGEV(p, loc = theta[1], scale = theta[2], shape = theta[3])
}

theta0 <- c(loc = mu, scale = sigma, shape = xi)
fval <- qGEV(p = p, loc = mu, scale = sigma, shape = xi,
             deriv = TRUE, hessian = TRUE)

Jnum <- jacobian(func = f, x = theta0)
Hnum <- hessian(func = f, x = theta0) 

test_that(desc = "Gradient of the GEV quantile",
            expect_lt(max(abs(attr(fval, "gradient") - Jnum)), 1e-6))

test_that(desc = "Hessian of the GEV quantile",
            expect_lt(max(abs(drop(attr(fval, "hessian")) - Hnum)),
                      1e-4))

## ==========================================================================
## Check the gradient and the Hessian of the distribution function
## ==========================================================================

n <- 20
mu <- rnorm(1, sd = 10)
sigma <- rexp(1)
xi <- rnorm(1, sd = 0.1)
x <- rGEV(n, loc = mu, scale = sigma, shape = xi) 

f <- function(theta) {
    pGEV(x, loc = theta[1], scale = theta[2], shape = theta[3])
}

theta0 <- c(loc = mu, scale = sigma, shape = xi)
fval <- pGEV(x, loc = mu, scale = sigma, shape = xi,
             deriv = TRUE)

Jnum <- jacobian(func = f, x = theta0)

test_that(desc = "Gradient of the GEV distribution function",
            expect_lt(max(abs(attr(fval, "gradient") - Jnum)), 1e-6))


