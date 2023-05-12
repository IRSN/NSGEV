##' @export
##' @method plot TVGEV
plot.TVGEV <- function(x, y, which = "c", ...) {
    Date <- loc <- NULL
    warning("This function is still in an early stage of\n",
            "developpement. Use the method 'plot' for bts\n",
            "or for predicted objects, and be patient!")

    
    if (which == "c") {
        indPar <- !x$isCst
        nPar <- sum(indPar)
        if (nPar == 0) stop("all GEV parameters are constant")
        theta <- data.frame(Date = x$data[ , x$date],
                            coef(x, type = "theta"))
        g <- ggplot(data = theta)
        g <- g + geom_line(mapping = aes(x = Date, y = loc))
    }
    g
}

##'
##' @title Create a \code{ggplot} for a \code{TVGEV} Object
##' 
##' @param object An object with class \code{"TVGEV"}.
##' 
##' @param type The type of result shown along with the series used.
##'     For now, only the timeseries of quantiles as computed by
##'     \code{\link{quantile.TVGEV}} can be added.
##' 
##' @param ... Not used yet
##'
##' @return A graphic object inheriting from \code{"ggplot"}.
##'
##' @method autoplot TVGEV
##' 
##' @note Mind that several methods of the class \code{"TVGEV"} such
##'     as \code{predict} \code{quantMax}, \code{coef} produce results
##'     that can be autoplotted.
##' 
autoplot.TVGEV <- function(object, type = c("quant"), ...) {

    Date <- Response <- value <- variable <- NULL
    
    type <- match.arg(type)

    df <- data.frame(Date = object$data[ , object$date],
                     Response = object$data[ , object$response])
    
    ymin <- min(df$Response, na.rm = TRUE)
    g <- ggplot(data = df)
    g <- g + geom_segment(mapping = aes(x = Date, xend = Date,
                                        y = Response, yend = ymin),
                          col = "black")

    if (type == "quant") {
        dfq <- as.data.frame(quantile(object))
        dfq <- tidyr::gather(dfq, key = type, value = value, -Date)
        g <- g + geom_line(data = dfq,
                           mapping = aes(x = Date, y = value, group = type,
                                         colour = type))
    }
    g <- g + ylab(object$response) + xlab("Date")
    g
}

