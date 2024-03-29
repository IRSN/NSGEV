##*****************************************************************************
##' Generalised Residuals for a \code{TVGEV} model.
##' 
##' @title Generalised Residuals for a \code{TVGEV} Model
##'
##' @aliases resid.TVGEV
##' 
##' @param object A \code{TVGEV} object.
##'
##' @param type The approximate distribution wanted. The choices
##'     \code{c("unif", "exp", "gumbel")} correspond to the standard
##'     uniform, the standard exponential and the standard Gumbel
##'     distributions. Partial matching is allowed.
##'
##' @param ... Not used yet.
##' 
##' @return A vector of generalised residuals which should
##' \emph{approximately} be independent and \emph{approximately}
##' follow the standard exponential or the uniform distribution,
##' depending on the value of \code{type}.
##' 
##' @note The upper 95\% quantile of the standard exponential is close
##'     to \eqn{3} which can be used to gauge "large residuals". Using
##'     \code{type = "gumbel"} seems better to diagnose abnormally
##'     small residuals as may result from abnormally small block
##'     maxima.
##'
##' @references Cox, D.R. and Snell E.J. (1968) "A General Definition
##' of Residuals".  \emph{JRSS Ser. B}, \bold{30}(2), pp. 248-275.
##'
##' Panagoulia, D. and Economou, P. and Caroni, C. (2014)
##' "Stationary and Nonstationary Generalized Extreme Value Modelling of
##'  Extreme Precipitation over a Mountainous Area under Climate Change".
##' \emph{Environmetrics} \bold{25}(1), pp. 29-43.
##' 
##' @importFrom nieve pGEV
##' @importFrom stats residuals resid
##' @method residuals TVGEV
##' @export
##' 
##' @examples
##' df <- within(TXMax_Dijon, Date <- as.Date(sprintf("%4d-01-01", Year)))
##' tv <- TVGEV(data = df, response = "TXMax", date = "Date",
##'             design = breaksX(date = Date, breaks = "1970-01-01", degree = 1),
##'             loc = ~ t1 + t1_1970)
##' e <- resid(tv)
##' plot(e)
##' ## ggplot alternative
##' autoplot(e)
##' ## plot the residual against the fitted location. Use 'as.numeric'
##' ## on the residuals to build a similar ggplot
##' mu <- tv$theta[ , "loc"]
##' plot(mu, e, type = "p", pch = 16, col = "darkcyan",
##'      main = "generalised residuals against 'loc'")
residuals.TVGEV <- function(object,
                            type = c("exp", "unif", "gumbel"),
                            ...) {

    type <- match.arg(type)
    Y <- object$data[ , object$response]
    theta <- psi2theta(model = object, psi = NULL)
    
    e <- nieve::pGEV(Y, loc = theta[, 1L], scale = theta[, 2L], 
                     shape = theta[, 3L], lower.tail = FALSE)

    lims95 <- c(0.025, 0.975)
    
    if (type == "exp") {
        e <- -log(e)
        lims95 <- -log(1 - lims95)
    } else if (type == "gumbel") {
        e <- log(-log(e))
        lims95 <- log(-log(1 - lims95))
    }
        
    names(e) <- rownames(theta)
    attr(e, "date") <- object$data[ , object$date]
    attr(e, "type") <- type
    attr(e, "lims95") <- lims95
    class(e) <- "resid.TVGEV"
    e
    
}

##*****************************************************************************
##' Plot the residuals of a \code{TVGEV} model object against the
##' date.
##'
##' @title Plot Residuals of a \code{TGVEV} Model
##'
##' @param x An object with class \code{"TVGEV"}.
##'
##' @param y Not used yet.
##'
##' @param ... Further arguments to be passed to \code{plot}.
##'
##' @return Nothing.
##'
##' @seealso \code{\link{residuals.TVGEV}}.
##'
##' @method plot resid.TVGEV
##' @export
##' 
plot.resid.TVGEV <- function(x, y = NULL, ...) {

    if (!missing(y) && !is.null(y)) {
        warning("'y' formal not used by this method")
    }
    type <- attr(x, "type")
    plot(attr(x, "date"), x, type = "o",
         pch = 16, col = "orangered",
         xlab = "date",
         ylab = sprintf("residual, type = \"%s\"", type),
         ...)
    abline(h = attr(x, "lims95"), col = "SpringGreen3")
    
}

##*****************************************************************************
##' Plot the residuals of a \code{TVGEV} model object against the
##' date.
##' 
##' @title Plot Residuals of a \code{TGVEV} Model
##'
##' @param object A \code{TVGEV} object.
##'
##' @param geom The "geometry" used to show the series of block
##'     maxima. The default \code{"line"} correspond to a polyline and
##'     \code{"segment"} correspond to a series of vertical segments
##'     a.k.a. "needles".
##' 
##' @param ... Not used yet.
##'
##' @return Nothing.
##'
##' @seealso \code{\link{residuals.TVGEV}}.
##'
##' @method autoplot resid.TVGEV
##' @export
##' 
autoplot.resid.TVGEV <- function(object,
                                 geom = c("line", "segment"),
                                 ...) {
    
    Date <- Resid <- NULL
    geom <- match.arg(geom)
    typeResid <- attr(object, "type")
    
    df <- data.frame(Date = attr(object, "date"), Resid = as.numeric(object))
    gg <- ggplot(data = df)
    
    if (geom == "line") {
        gg <- gg + geom_line(mapping = aes(x = Date, y = Resid),
                             colour = "orangered")
    } else {
        gg <- gg + geom_segment(mapping = aes(x = Date, xend = Date,
                                              y = 0, yend = Resid),
                                colour = "orangered")
    }
    gg <- gg + geom_point(mapping = aes(x = Date, y = Resid),
                          colour = "orangered", shape = 16)

    gg <- gg + geom_hline(data = data.frame(lims = attr(object, "lims95")),
                          mapping = aes(yintercept = lims),
                          col = "SpringGreen3", alpha = 0.7)
    gg <- gg + xlab("date") + ylab(sprintf("residual type = \"%s\"", typeResid)) +
        ggtitle("Generalized residuals with 95% limits")
    gg
    
}
