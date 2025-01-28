##' This function provides some indications on the influence of a
##' specific observation used in the fit of a \code{TVGEV} object. For
##' several candidate values of the chosen observation, a \code{TVGEV}
##' is re-fitted with the observed value replaced by the candidate
##' value, and the quantity of interest as defined by \code{what} is
##' computed.
##' 
##' By suitably defining the function \code{what}, one can trace e.g.,
##' the estimate of a quantile with given probability, the
##' corresponding upper confidence limit, the upper end-point of the
##' distribution, ... and more.
##' 
##' @title Diagnostic of Influence for \code{TVGEV} Objects
##' 
##' @param model A \code{TVGEV} object. 
##'
##' @param what A function with one argument \code{object} defining or
##'     extracting the (numeric) quantity for which the influence is
##'     to be analysed. So, \code{object} must be a \code{TVGEV}
##'     object. See \bold{Examples}. This function can return numeric
##'     vector in which case several influence curves will be
##'     computed. It is then wise to make the function return a vector
##'     with named elements.
##'
##' @param which Identifies the observation which influence will be
##'     analysed. Can be: an integer giving the observation, the
##'     character \code{"min"}, the character{"max"}, a character that
##'     can be coerced to a date, or an object with class
##'     \code{"Date"} that identifies the observation.  This argument
##'     must define an observation existing in the data frame used to
##'     create the object. So it can not be used e.g. to define a
##'     future date.
##'
##' @param how Character specifying the type of (finite-sample)
##'     influence curve to use. For now only the only possibility is
##'     \code{"replace"}, see Section 2.1 of the book by Hampel et al
##'     (1996).  The chosen observation is replaced by a value \eqn{x}
##'     and the quantity of interest is computed for each value of
##'     \eqn{x}.
##'
##' @param ... Not used.
##'     
##' @return A named list of vectors. The element \code{obsValue}
##'     contains the candidate values for the chosen observation and
##'     \code{statValue} contains the corresponding values of the
##'     specific statistic defined with \code{what}.
##'
##' @note For \code{TVGEV} models with a positive GEV shape, the
##'     smallest observation may have a strong influence on the
##'     result, even if the model is stationary.
##'
##' @importFrom stats influence
##' @export
##' @method influence TVGEV
##' 
##' @references
##'
##'   Hampel, F.R., Ronchetti, E.M., Rousseeuw, P.J. and Stahel
##'   W.A. (1996) \emph{Robust Statistics: The Approach Based on
##'   Influence Functions}. Wiley.
##' 
##' @examples
##' example(TVGEV)
##' influ <- influence(res2, which = "min")
##' autoplot(influ)
##' influ <- influence(res2, which = as.Date("2015-01-01"))
##' autoplot(influ)
##' RL_2050 <- function(model) {
##'    c("RL_2050(100)" = quantile(model, prob = 0.99, date = "2050-01-01")[1])
##' }
##' influence(res2, what = RL_2050) |> autoplot()
##' influence(res2, what = RL_2050, which = "2003-01-01") |> autoplot()
##' RLs_2050 <- function(model) {
##'    res <- c(quantile(model, prob = 0.99, date = "2050-01-01")[1],
##'             quantile(model, prob = 0.999, date = "2050-01-01")[1])
##'    names(res) <- paste0("RL", c(100, 1000))
##'    res
##' }
##' influ <- influence(res2, what = RLs_2050, which = "2003-01-01")
##' autoplot(influ)
influence.TVGEV <- function(model,
                            what = function(model) coef(model)["xi_0"],
                            which = "min",
                            how = "replace",
                            ...) {
    
    y <- model$data[ , model$response]
    Date <- model$data[ , model$date]
    n <- model$n
    ind <- (1:n)[!is.na(y)]
    if (is.numeric(which)) {
        if (abs(which - round(which)) > 0.01 || which < 1 || which > n) {
            stop("When 'which' is numeric, it must be an integer between ",
                 "1 and `nobs(model)`")
        }
    } else if (is.character(which)) {
        if (which == "min") {
            ii <- which.min(y[ind])
            which <- ind[ii]
        } else if (which == "max") {
            ii <- which.max(y[ind])
            which <- ind[ii]
        } else { 
            which <- as.Date(which)
            which  <- (1:n)[Date == which]
        }
    } else if (inherits(which, "Date")) {       
        which <- (1:n)[Date == which]
    }

    rr <- range(y, na.rm = TRUE)
    rr <- rr + c(-0.1, 0.1) * diff(rr)

    yGrid <- seq(from = rr[1], to = rr[2], length = 30)
    yGrid <- c(yGrid, y[which])
    res0 <- do.call(what, list(model = model))
    nres <- length(res0)
    nms <- names(res0)
    res <- matrix(NA_real_, nrow = length(yGrid), ncol = nres)
    names(res) <- nms
    
    for (iy in seq_along(yGrid)) {
        .df <- model$data
        .df[which, model$response] <- yGrid[iy]
        calli <- model$call
        calli[["data"]] <- .df
        tr <- try(fiti <- eval(calli))
        resi <- do.call(what, list(model = fiti))
        res[iy, ] <- resi
    }

    if (nres == 1) res <- as.vector(res)
    
    res <- list(obsValue = yGrid,
                statValue = res,
                which = which,
                what = what,
                date = Date,
                y = y,
                names = nms)
    
    class(res) <- "influence.TVGEV"
    res
    
}

##' @method print influence.TVGEV
##' @export
##' 
print.influence.TVGEV <- function(x, ...) {
    cat("Finite-sample influence function\n")
    cat(sprintf("Observation number: %d (%s)\n",
                   x$which, x$date[x$which]))
    cat(sprintf("what:               %s\n", paste(x$names, collapse = ", ")))
}

##'
##' @title Plot Method for \code{influence.TVGEV} Object 
##'
##' @param x A \code{TVGEV} object.
##'
##' @param y Not used.
##'
##' @param ... Not used.
##'
##' @return Nothing.
##'
##' @note This method \code{autoplot}
##'
##' @seealso \code{\link{autoplot.TVGEV}}.
##'
##' @importFrom grDevices palette adjustcolor
##' @importFrom graphics matplot grid rug legend
##' @export
##' @method plot influence.TVGEV
##' 
plot.influence.TVGEV <- function(x, y = NULL, ...) {
    Pch <- c(16, 18, 21, 24, 22)
    nc <- ifelse(is.matrix(sv <- x$statValue), ncol(sv), 1)
    warning("Consider using the `autoplot` method instead")
    grDevices::palette("Tableau")
    matplot(x$obsValue, x$statValue,
            pch = Pch[1:nc],
            main = sprintf("Influence of obs. number %d, (%s)",
                           x$which, format(x$date[x$which])),
            xlab = "obs. candidate value", ylab = "stat. value")
    grid()
    abline(v = x$y[x$which],
           col = grDevices::adjustcolor( "gray", alpha.f = 0.7))
    rug(x$y)
    legend("topleft",
           pch = Pch[1:nc],
           col = 1:nc,
           legend = x$names)
    grDevices::palette("default")
}


##' The influence curve is obtained by replacing the chosen
##' observation by a candidate value and then re-fitting the same
##' model and computing the statistic of interest: coefficient
##' estimate, quantile, ... The rug shows the \eqn{y_b} and the
##' vertical line shows the observation chosen for the analysis.
##' 
##' @title Create a \code{ggplot} for a \code{influence.TVGEV} Object 
##' 
##' @param object An object with class \code{"influence.TVGEV"} as
##'     created bye the \code{influnce} method for the \code{"TVGEV"}
##'     class.
##' 
##' @param ... Not used yet.
##' 
##' @return An object inheriting grom \code{"ggplot"} showing the
##'     (finite-sample) influence function for the observation and the
##'     statistic that have been chosen.
##'
##' @export
##' @method autoplot influence.TVGEV
##'
##' @examples
##' library(ismev)
##' data(venice)
##' df <- data.frame(Date = as.Date(paste0(venice[ , "Year"], "-01-01")),
##'                  Sealevel = venice[ , "r1"] / 100)
##' fit0 <- TVGEV(data = df, date = "Date", response = "Sealevel")
##' autoplot(fit0)
##' RL_2050 <- function(model) {
##'     c("RL_2050(100)" = quantile(model, prob = 0.99, date = "2050-01-01")[1])
##' }
##' autoplot(fit0)
##' influence(fit0, what = RL_2050) |> autoplot()
##' ## fit with a linear time trend
##' fit1 <- TVGEV(data = df, date = "Date", response = "Sealevel",
##'               design = polynomX(Date, degree = 1), loc = ~t1)
##' autoplot(fit1)
##' summary(fit1)
##' influence(fit1) |> autoplot()
##' influence(fit1, what = RL_2050) |> autoplot()
##' influence(fit1, which = "max", what = RL_2050) |> autoplot()
##' ## Influence curve for the estimated slope
##' slope <- function(model) {
##'     c("slope\n(cm/100 yr)" = 100 * unname(coef(model)["mu_t1"]))
##' }
##' influence(fit1, which = "max", what = slope) |> autoplot()
##' influence(fit1, which = "min", what = slope) |> autoplot()
autoplot.influence.TVGEV <- function(object, ...) {

    obsValue <- statValue <- name <- y <- NULL
    
    if (is.matrix(m <- object$statValue) && ncol(m) > 1) { 
        df <- list()
        for (j in 1:ncol(object$statValue)) {
            df[[j]] <- data.frame(obsValue = object$obsValue,
                                  statValue = object$statValue[ , j],
                                  name = object$names[j])
        }
        df <- data.table::rbindlist(df)
    } else {
        df <- data.frame(obsValue = object$obsValue,
                         statValue = object$statValue,
                         name = object$names)
    }
    g <- ggplot(data = df) +
        geom_point(mapping = aes(x = obsValue,
                                 y = statValue,
                                 group = name,
                                 colour = name,
                                 shape = name), size = 1.5) +
        xlab("obs. value") + ylab("stat. value") +
        geom_vline(xintercept = object$y[object$which],
                   colour = adjustcolor( "gray", alpha.f = 0.9)) +
        ggtitle(sprintf("Influence curve of obs. number %d, (%s)",
                        object$which, format(object$date[object$which]))) + 
        geom_rug(data = data.frame(y = object$y), mapping = aes(x = y)) +
        labs(colour = "stat.", shape = "stat.") +
        scale_colour_manual(values = c("orangered", "springgreen3", "steelblue3"))
    g
    
}
