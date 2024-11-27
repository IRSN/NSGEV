## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## GOAL: Test the gradient and the Hessian as computed by the
## `qMax.TVGEV` function. This function computes the quantile of the
## maximum 'M' on a given period of time from a `TVGEV` model object.
##
## ***************************************************************************

library(numDeriv)
library(testthat)
context("qMax.TVGEV")

HESSIAN <- FALSE

library(NSGEV)
df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))

## =============================================================================
## - fit a TVGEV model. Only the location parameter is TV.
## 
## - Define several scenarios for the "design life" period on which
##   the max is to be computed and the probability for the quantile.
##
## - Define several probabilities (non-exceedance) for which the
##   quantile of the maximum will be computed.
##   =============================================================================

tv  <- TVGEV(data = df, response = "TXMax", date = "Date",
             design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
             loc = ~ t1_1970)
newDates <- list("one year" = as.Date(sprintf("%4d-01-01", 2025)),
                 "ten years" = as.Date(sprintf("%4d-01-01", 2025:2034)),
                 "thirty years" = as.Date(sprintf("%4d-01-01", 2025:2054)))
p <- c(0.75, 0.8, 0.9, 0.95, 0.97, 0.99, 0.995)

## =============================================================================
## For each scenario and each probability, compute the quantiles with
## its gradient and its Hessian. Compare the value of these with the
## corresponding derivatives computed by numeric differentiation.
##
## For now the Hessian part of the check is not run, because the test
## fails for high probabilities and long periods. There is some
## evidence that the symbolic Hessian is exact, yet that the numeric
## computation has a poor precision. Yet the computation of the
## symbolic Hessian could also be subject to problems of numeric
## precision as well.
## =============================================================================

for (iDate in seq_along(newDates)) {
    for (ip in seq_along(p)) {
        
        pM <- qMax.TVGEV(object = tv, p = p[ip],
                         date = newDates[[iDate]], deriv = TRUE, hessian = HESSIAN)
        
        myFun <- function(psi) {
            qMax.TVGEV(object = tv, p = p[ip], date = newDates[[iDate]],
                       psi = psi)
        }

        psi0 <- coef(tv)

        gradNum <- numDeriv::grad(func = myFun, x = psi0)
        grad <- drop(attr(pM, "gradient"))
        
        cond <- testNumDeriv(grad, gradNum, type = "gradient")
        
        if (!cond) {
            cat(sprintf("\nGradient of qMax %s, p = %5.4f\n",
                        names(newDates)[iDate], p[ip]))
            cat("Provided and numeric gradients\n")  
            print(cbind("provided" = round(grad, digits = 4),
                  "numeric" = round(gradNum, digits = 4))) 
        }
        test_that(desc = sprintf("Gradient of qMax %s, p = %5.4f",
                                 names(newDates)[iDate], p[ip]),
                  expect_true(cond))

        if (HESSIAN) {
            hess <- drop(attr(pM, "hessian"))
            hessNum <-  numDeriv::hessian(func = myFun, x = psi0)
            cond <- testNumDeriv(hess, hessNum, type = "hessian")
            
            if (!cond) {
                cat(sprintf("\nHessian of qMax %s, p = %5.4f\n",
                            names(newDates)[iDate], p[ip]))
                cat("Provided and numeric Hessians\n")
                dimnames(hessNum) <- dimnames(hess)
                print(round(hess, dig = 4))
                print(round(hessNum, dig = 4))
                trace <- 1
                cat(sprintf("Re-run with 'trace = %d'\n", trace))
                pM <- qMax.TVGEV(object = tv,
                                 p = p[ip],
                                 date = newDates[[iDate]],
                                 deriv = TRUE, hessian = TRUE,
                                 trace = trace)
            } else {
                cat(sprintf("\nHessian of qMax %s, p = %5.4f: OK!\n",
                            names(newDates)[iDate], p[ip]))
            }
            test_that(desc = sprintf("Hessian of qMax %s, p = %5.4f",
                                     names(newDates)[iDate], p[ip]),
                      expect_true(cond))
        }
    }
}
    
