## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## GOAL: Test the gradient as computed by the `dMax.TVGEV`
## function. Remind that this function computes the density function
## of the maximum 'M' on a given period of time from a `TVGEV` model
## object.
##
## NOTE: For now the Hessian is not computed by the function. The
## derivative w.r.t. 'x' as computed by the function will not be
## checked as this is a very simple computation. Also the derivative
## w.r.t. 'x' is used in the Hessian of the quantile `qMax.TVGEV`,
## hence no sucessful test for the Hessian could be obtained is the
## result provided by `dMax.TVGEV` was wrong.
##
## ***************************************************************************

library(numDeriv)
library(testthat)
context("dMax.TVGEV")

HESSIAN <- TRUE
TRACE <- FALSE

library(NSGEV)
df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))

## =============================================================================
## - fit a TVGEV model. Only the location parameter is TV.
## 
## - Define several scenarios for the "design life" period on which
##   the max is to be computed and the probability for the quantile.
##
## - Define several quantiles for which the distribution function of
##   the maximum will be computed.
## =============================================================================

tv  <- TVGEV(data = df, response = "TXMax", date = "Date",
             design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
             loc = ~ t1_1970)
newDates <- list("one year" = as.Date(sprintf("%4d-01-01", 2025)),
                 "ten years" = as.Date(sprintf("%4d-01-01", 2025:2034)),
                 "thirty years" = as.Date(sprintf("%4d-01-01", 2025:2054)))
x <- c(40.0, 42.0, 44.0, 46.0, 48.0, 50.0)

## =============================================================================
## For each scenario and each quantile, compute the value of the
## density with its gradient. Compare the value of this with the
## corresponding derivatives computed by numeric differentiation.
## =============================================================================

for (iDate in seq_along(newDates)) {
    for (ix in seq_along(x)) {
        
        dM <- dMax.TVGEV(object = tv, x = x[ix],
                         date = newDates[[iDate]], deriv = TRUE, hessian = HESSIAN)
        
        myFun <- function(psi) {
            dMax.TVGEV(object = tv, x = x[ix], date = newDates[[iDate]],
                       psi = psi)
        }
        
        psi0 <- coef(tv)
        
        gradNum <- numDeriv::grad(func = myFun, x = psi0)
        grad <- drop(attr(dM, "gradient"))

        if (TRACE) print(cbind(provided = grad, num = gradNum))
        
        cond <- testNumDeriv(grad, gradNum, type = "gradient")
        
        if (!cond) {
            cat(sprintf("\nGradient of qMax %s, x = %5.4f\n",
                        names(newDates)[iDate], x[ix]))
            cat("Provided and numeric gradients\n")  
            print(cbind("provided" = round(grad, digits = 4),
                  "numeric" = round(gradNum, digits = 4))) 
        }
        test_that(desc = sprintf("Gradient of qMax %s, x = %5.4f",
                                 names(newDates)[iDate], x[ix]),
                  expect_true(cond))
    }
}
    
