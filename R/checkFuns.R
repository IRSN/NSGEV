## ============================================================================
##' Check the \code{confint method}.
##'
##' Confidence intervals on the parameters of \code{object} are
##' computed using the \code{confint} method for the class
##' \code{"TVGEV"}. The confidence interval for the shape parameter is
##' compared to that given by the \code{ismev::prof.xi} function.
##'
##' @title Check the confint Method for \code{TVGEV} Objects
##'
##' @param object An object with class \code{"TVGEV"} representing a
##' stationary model.
##'
##' @param ref An object with class \code{"gev.fit"} created by the
##' function \code{gev.fit} of the \pkg{ismev} package. 
##'
##' @param level Confidence level.
##'
##' @param which Not used yet.
##'
##' @param method See \code{\link{confint.TVGEV}}. For now, only
##' \code{"proflik"} is possible.
##'
##' @return Nothing is returned: on the plot built by
##' \code{ismev::prof.xi}, two thick vertical lines are added, showing
##' the confidence limits found by \code{confint.TVGEV} with
##' \code{method = "proflik"}. Each vertical line should cut
##' the profiled log-likelihood curve at one of its two intersections
##' with the horizontal line.
##'
##' @seealso \code{\link{confint.TVGEV}}.
##'
##' @export
##' 
##' @examples
##' data(portpirie)
##' df <- portpirie
##' df <- within(df, date <- as.Date(sprintf("%4d-01-01", Year)))
##' fitNSGEV <- try(TVGEV(data = df,
##'                        response = "SeaLevel", date = "date"))
##' fitismev <- gev.fit(xdat = df$SeaLev)
##' check_confint(fitNSGEV, ref = fitismev)
##'
check_confint <- function(object,
                          ref,
                          level = 0.95,
                          which = "shape",
                          method = "proflik") {

    if (missing(ref) || !inherits(ref, "gev.fit")) {
        stop("`ref' must be an object with class \"gev.fit\"",
             "created by using the 'ismev::fevd' function")
    }
    
    ciTVGEV <- confint(object, level = level, method = "proflik")
    
    rbind(coef(object), ref$mle)
    gev.profxi(ref, conf = level, xlow = -0.5, xup = 0.5)
    
    abline(v = ciTVGEV["xi_0", , 1], col = "SpringGreen3", lwd = 2)
    if (FALSE) {
        abline(v = c(-0.241, 0.142), col = "orangered", lwd = 2)
        title(main = "profile likelihood for shape xi. Typo in BOOK.")
        legend("bottomright",
               lty = "solid", lwd = 2,
               col = c("orangered", "SpringGreen3"),
               legend = c("ismev BOOK CI", "NSGEV CI"))
    }
}

## ============================================================================
##' Check the \code{predict method}.
##'
##' The comparison with the \pkg{ismev} package is only completed for
##' stationary models because it is difficult otherwise to get
##' profile-likelihood confidence intevals on Return Levels.
##'
##' @title Check the Predict Method for \code{TVGEV} Objects
##'
##' @param object A \code{TVGEV} package.
##'
##' @param ref An object with class \code{"gev.fit"} created by the
##' function \code{gev.fit} of the \pkg{ismev} package. 
##'
##' @param level Confidence level.
##' 
##' @param which Not used yet.
##'
##' @param newdate See \code{\link{predict.TVGEV}}.
##'
##' @param confintMethod See \code{\link{predict.TVGEV}}.
##'
##' @return Nothing. For each return period computed by
##' \code{predict.TVGEV}, the plot built by \code{ismev::gev.prof} is
##' shown. Two vertical lines are added, showing
##' the confidence limits found by \code{predict.TVGEV} with
##' \code{confintMethod = "proflik"}. Each vertical line should cut
##' the profiled log-likelihood curve at one of its two intersections
##' with the horizontal line.
##'
##' @note \code{ismev::gev.prof} does not return the confidence limits
##' which must be evaluated by eye from the displayed graph.
##'
##' @seealso \code{\link{predict.TVGEV}}.
##'
##' @importFrom graphics text par
##' 
##' @export
##' 
##' @examples
##' data(portpirie)
##' df <- portpirie
##' df <- within(df, date <- as.Date(sprintf("%4d-01-01", Year)))
##' fitNSGEV <- try(TVGEV(data = df,
##'                        response = "SeaLevel", date = "date"))
##' fitismev <- gev.fit(xdat = df$SeaLev)
##' check_predict(fitNSGEV, ref = fitismev)
##'
check_predict <- function(object,
                          ref,
                          level = 0.95,
                          which = "confint",
                          newdate = as.Date("2020-02-01"),
                          confintMethod = "proflik") {
    
    Period <- NULL ## avoid WARNING in package check
    
    if (missing(ref) || !inherits(ref, "gev.fit")) {
        stop("`ref' must be an object with class \"gev.fit\"",
             "created by using the 'ismev::fevd' function")
    }
    
    predTVGEV <- predict(object,
                         newdate = newdate,
                         level = level,
                         confintMethod = confintMethod)
    
    if (confintMethod == "proflik") {
        
        opar <- par(mfrow = c(ceiling(nrow(predTVGEV) /3) , 3),
                    mar = c(2, 5, 1, 1), oma = c(1, 1, 1, 1))
        
        for (m in unique(predTVGEV$Period)) {
            LU <- subset(predTVGEV, Period == m)[ , c("L", "U")]
            LU <- as.numeric(LU)
            LUmod <- LU + c(-1, 1) * diff(LU) / 20
            gev.prof(ref, m = m,
                     conf = level,
                     xlow = LUmod[1], xup = LUmod[2])
            abline(v = LU,
                   col = "SpringGreen3", lwd = 2)
            pu <- par()$usr
            text(x = (pu[1] + pu[2]) / 2, y =  (pu[3] + pu[4]) / 2,
                 labels = sprintf("Period = %d", m),
                 col = "SpringGreen4", cex = 1.4)
        }
        
        par(opar)
        
    }
    
}

## ============================================================================
##' Check the 'predictUncond' function.
##'
##' The unconditional prediction is computed from the object given
##' using \code{predictUncond} this leads to a table, say
##' \code{pu}. Then \code{nsim} paths are simulated from the model;
##' The simulation is for the period of time begining at
##' \code{newdateFrom} and having the length suitable to cover all the
##' predicted periods found in \code{pu}. For each row of the
##' prediction table \code{pu} -by default, corresponding to a period 5, 10, 20,
##' ...- we count the number of exceedances over the quantile given in
##' the \code{'Quant'} column, for each sample. By definition of the
##' Non-Stationary Return Level, the random number of exceedances has
##' an expectation of one. The function returns the average numbers of
##' exceedances over the \code{nsim} paths which therefore should be
##' close to one at least if \code{nsim} is large.
##' 
##' @title Check the \code{\link{predictUncond}} Function
##'
##' @param object An object with class \code{"TVGEV"}
##'
##' @param newdateFrom Argument passed to \code{\link{predictUncond}}.
##'
##' @param nsim Number os simulations to use in in the check.
##'
##' @param plot Logical. If \code{TRUE}, some simulated paths will be
##' shown.
##' 
##' @param trace Integer level of verbosity.
##' 
##' @param ... Other arguments to be passed to \code{\link{predictUncond}}.
##'
##' @return A vector of average number of exceedances corresponding to
##' each of the periods which have been predicted. The elemnts of this
##' vector should be close to one, see \bold{Details}.
##'
##' @seealso \code{\link{predictUncond}}.
##'
##' @export
##'
##' @examples
##' example(TVGEV)
##' set.seed(1357)
##' check_predictUncond(res2)
##' check_predictUncond(res2, newdateFrom = "1900-01-01", period = 200)
##' 
check_predictUncond <- function(object,
                                newdateFrom = "2020-01-01",
                                nsim = 5000,
                                plot = TRUE,
                                trace = 0,
                                ...) {
    
    pu <- predictUncond(object, newdateFrom = newdateFrom, ...)
    if (trace) {
        cat("Unconditional predictions\n")
        print(pu)
    }
    
    period <- pu$Period
    perPred <- seq(from = as.Date(newdateFrom), by = "year",
                  length.out = ceiling(max(period)))
    ysim <- simulate(object, nsim = 6000, newdate = perPred)
    
    if (plot) {
        nsimMax <- 100
        if (nsim > nsimMax) {
            ysim2 <- ysim[ , 1:nsimMax]
            attr <- attributes(ysim)
            attr[["dim"]] <- c(attr[["dim"]][1], nsimMax)
            attr[["dimnames"]][[2]] <- attr[["dimnames"]][[2]][1:nsimMax] 
            attributes(ysim2) <- attr
            plot(ysim2, main = sprintf("First %d simulations", nsimMax))
        } else {
            plot(ysim)
        }
    }
    
    period <- pu$Period
    rho <- pu$Quant
    rp <- round(period)
    if (any((period - rp) > 1e-4)) {
        stop("'period' must contain integer values")
    } else period <- rp

    nExceed <- matrix(NA, nrow = length(rho), ncol = ncol(ysim))
    for (i in seq_along(rho)) {
        nExceed[i, ] <- apply(ysim[1:period[i], ], 2,
                              function(x) { sum(x > rho[i]) })
    }
    
    check <- apply(nExceed, 1, mean)
    names(check) <- period
    check
    
}

