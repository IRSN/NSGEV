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
##' @param geom Character giving the "geometry" to be used for the
##'     layer corresponding to the (timeseries) data.
##' 
##' @param type The type of result shown along with the series used.
##'     For now, only the timeseries of quantiles as computed by
##'     \code{\link{quantile.TVGEV}} can be added.
##' 
##' @param ... Not used yet
##'
##' @return A graphic object inheriting from \code{"ggplot"}.
##' 
##' @note Mind that several methods of the class \code{"TVGEV"} such
##'     as \code{predict} \code{quantMax}, \code{coef} produce results
##'     that can be autoplotted.
##'
##' @method autoplot TVGEV
##' @export
##' 
autoplot.TVGEV <- function(object, geom = c("line", "segment"),
                           type = c("quant"), ...) {
    
    Date <- Response <- value <- variable <- NULL
    dots <- match.call(expand.dots = FALSE)$...
    geom <- match.arg(geom)
    type <- match.arg(type)
    colour <- "orangered"
    
    if (!is.na(m <- match("colour", names(dots)))) colour <- dots$colour
        
    df <- data.frame(Date = object$data[ , object$date],
                     Response = object$data[ , object$response])
    
    ymin <- min(df$Response, na.rm = TRUE) -
        diff(range(df$Response, na.rm = TRUE)) / 10
    g <- ggplot(data = df)
    
    if (geom == "segment") {
        g <- g + geom_segment(mapping = aes(x = Date, xend = Date,
                                            y = ymin, yend = Response),
                              col = colour)
    } else {
        g <- g + geom_line(mapping = aes(x = Date, y = Response),
                           col = colour)
    }
    
    g <- g + geom_point(mapping = aes(x = Date, y = Response),
                        shape = 16, col = colour)
    
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

