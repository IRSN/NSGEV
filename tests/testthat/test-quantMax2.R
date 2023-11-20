# ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
##
## GOAL: Test that the quantile of the maximum is correctly computed
## by comparing the result to that of an alternative computation based
## on a large number of simulations. Using a large number of simulated
## maxima M_i we can compute the related Empirical Cumulated
## Distribution Function (ecdf) say HatF(x). Then for each computed
## quantile 'q(p)', the value HatF{q(p)} should be close enough to the
## corresponding probability 'p'.
##
## ***************************************************************************

library(NSGEV)
library(testthat)
context("TVGEV: Quantile of the maximum: compare with simulations")

df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
tv <- TVGEV(data = df, response = "TXMax", date = "Date",
            design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
            loc = ~ t1_1970)

newDate <- as.Date(sprintf("%4d-01-01", 2025:2055))
qM <- quantMax(tv, date = newDate, level = 0.95)

nSim <- 1e5
sim <- simulate(tv, nsim = nSim, newdate = newDate)
M <- apply(sim, 2, max)
prob <- qM$Prob
qMSim <- quantile(M, prob = prob)

delta <- ecdf(M)(qM$Quant) - prob
s <- 4.0 * sqrt(prob * (1 - prob)) / sqrt(nSim)

test_that(desc = "The ecdf of the quantiles of Max should match the prob",
          expect_lt(max(abs(delta / s)), 1.0))
