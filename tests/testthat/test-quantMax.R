## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## GOAL: Test that the lower and upper confidence bounds for the
## quantile of the maximum are OK, be the 'date' argument used or
## not. Since the computation involves using the implicit function
## theorem we check that that we get nearly the same results is a
## numeric gradient is used on the quantile function, here assumed to
## be exact.
##
## ***************************************************************************

library(NSGEV)
library(numDeriv)
library(testthat)
context("TVGEV: Quantile of the maximum")

df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
tv <- TVGEV(data = df, response = "TXMax", date = "Date",
            design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
            loc = ~ t1 + t1_1970)


date2 <- as.Date(sprintf("%4d-01-01", 2025:2055))
qM1 <- quantMax(tv, level = 0.95)
qM2 <- quantMax(tv, date = date2, level = 0.95)

psiHat <- coef(tv)
Cov <- vcov(tv)

quantMaxPsi1 <- function(psi, prob) {
    qf <- quantMaxFun(object = tv, date = NULL, psi = psi)
    qf(prob)
}

quantMaxPsi2 <- function(psi, prob) {
    qf <- quantMaxFun(object = tv, date = date2, psi = psi)
    qf(prob)
}

q1 <- quantMaxPsi1(psiHat, qM1$Prob)
e <- q1- qM1$Quant 
test_that(desc = "Quantile Value. Case with 'date' unused",
          expect_lt(max(abs(e)), 5e-3))

q2 <- quantMaxPsi2(psiHat, qM2$Prob)
e <- q2- qM2$Quant 
test_that(desc = "Quantile Value. Case with 'date' used",
          expect_lt(max(abs(e)), 5e-3))

Var1 <- rep(0.0, nrow(qM1))
for (i in 1:nrow(qM1)) {
    J1 <- numDeriv::grad(func = quantMaxPsi1,
                         x = psiHat, prob = qM1$Prob[i])
    Var1[i] <- t(J1) %*% Cov %*% J1 
}

e <- qM1$Quant +  qnorm(0.975) * sqrt(Var1) - qM1$U 
test_that(desc = "Quantile Gradient w.r.t. 'psi'. Case with 'date' unused",
          expect_lt(max(abs(e)), 5e-3))

Var2 <- rep(0.0, nrow(qM2))
for (i in 1:nrow(qM2)) {
    J2 <- numDeriv::grad(func = quantMaxPsi2,
                         x = psiHat, prob = qM2$Prob[i])
    Var2[i] <- t(J2) %*% Cov %*% J2 
}

e <- qM2$Quant +  qnorm(0.975) * sqrt(Var2) - qM2$U 
test_that(desc = "Quantile Gradient w.r.t.t 'psi'. Case with 'date' used",
          expect_lt(max(abs(e)), 5e-3))
