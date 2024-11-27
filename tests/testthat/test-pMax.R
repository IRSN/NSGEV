## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## GOAL: Test the gradient and the Hessian as computed by the
## `pMax.TVGEV` function. This function computes the distribution
## function of the maximum 'M' on a given period of time from a
## `TVGEV` model object.
##
## ***************************************************************************

library(numDeriv)
library(testthat)
context("pMax.TVGEV")

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
q <- c(40.0, 42.0, 44.0, 46.0, 48.0)

## =============================================================================
## For each scenario and each probability, compute the quantiles with
## its gradient and its Hessian. Compare the value of these with the
## corresponding derivatives computed by numeric differentiation.
##
## The Hessian part is the tickiest!
## =============================================================================

for (iDate in seq_along(newDates)) {
    for (iq in seq_along(q)) {
        
        pM <- pMax.TVGEV(object = tv, q = q[iq],
                         date = newDates[[iDate]], deriv = TRUE, hessian = HESSIAN)
        
        myFun <- function(psi) {
            pMax.TVGEV(object = tv, q = q[iq], date = newDates[[iDate]],
                       psi = psi)
        }
        
        psi0 <- coef(tv)
        
        gradNum <- numDeriv::grad(func = myFun, x = psi0)
        grad <- drop(attr(pM, "gradient"))
        
        cond <- testNumDeriv(grad, gradNum, type = "gradient")
        
        if (!cond) {
            cat(sprintf("\nGradient of qMax %s, p = %5.4f\n",
                        names(newDates)[iDate], q[iq]))
            cat("Provided and numeric gradients\n")  
            print(cbind("provided" = round(grad, digits = 4),
                  "numeric" = round(gradNum, digits = 4))) 
        }
        test_that(desc = sprintf("Gradient of qMax %s, q = %5.4f",
                                 names(newDates)[iDate], q[iq]),
                  expect_true(cond))

        if (HESSIAN && TRACE) {
            hess <- drop(attr(pM, "hessian"))
            hessNum <-  numDeriv::hessian(func = myFun, x = psi0)
            cond <- testNumDeriv(hess, hessNum, type = "hessian")
            
            if (!cond) {
                cat(sprintf("\nHessian of qMax %s, q = %5.4f\n",
                            names(newDates)[iDate], q[iq]))
                cat("Provided and numeric Hessians\n")
                dimnames(hessNum) <- dimnames(hess)
                print(round(hess, dig = 4))
                print(round(hessNum, dig = 4))
                trace <- 1
                cat(sprintf("Re-run with 'trace = %d'\n", trace))
                pM <- qMax.TVGEV(object = tv,
                                 p = q[iq],
                                 date = newDates[[iDate]],
                                 deriv = TRUE, hessian = TRUE,
                                 trace = trace)
            } else {
                cat(sprintf("\nHessian of qMax %s, q = %5.4f: OK!\n",
                            names(newDates)[iDate], q[iq]))
            }
            test_that(desc = sprintf("Hessian of qMax %s, q = %5.4f",
                                     names(newDates)[iDate], q[iq]),
                      expect_true(cond))
        }
    }
}
    
