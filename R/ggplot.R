## ****************************************************************************
##' Plot a block time-series using ggplot.
##'
##' @title Plot a Block Time-Series using ggplot
##'
##' @param object An object with class \code{"bts"}.
##'
##' @param col1 A colour to be used when the number of component
##' series in \code{object} is > 12.
##'
##' @param facets Logical. If true, each component series will be
##' shown in a facet.
##' 
##' @param ... Not used
##'
##' @return An object of class \code{"ggplot"}.
##'
##' @import ggplot2
##' @importFrom forecast autoplot autolayer
##' @method autoplot bts
##' @export
##' 
autoplot.bts <- function(object, col1 = "darkgray", facets = FALSE, ...) {
    
    value <- variable <- type <- Date <- NULL
    nc <- ncol(object)
    ylab <- attr(object, "label")
    if (is.null(ylab)) ylab <- ""
    col <- "orangered"
    alpha <- opacity(nc)
    

    if (((nc == 1L) || (nc > 12L)) && facets) {
        warning("Value of 'facet' ignored because the number of ",
                "columns is 1 or is > 12")
    }
    
    g <- ggplot()
    
    if (nc == 1L) {
        df <- data.frame(Date = as.Date(attr(object, "date")),
                         value = unclass(object)[ , 1L])
        g <- g + geom_line(data = df,
                           mapping = aes(x = Date, y = value),
                           alpha = alpha, colour = col)
        g <- g + labs(colour = ylab, x = "", y = ylab)
    } else {
        
        if (FALSE) {
            df <- data.frame(date = as.Date(attr(object, "date")), object)
            df <- reshape2::melt(df, id.vars = "date")
        } else {
            df <- as.data.frame(object, longFormat = TRUE)
            ## df <- reshape2::melt(df, id.vars = "Date")
        }
        
        clab <- attr(object, "collabels")
        if (!is.null(clab)) levels(df$type) <- clab
        
        if (nc <= 12L) {
            g <- g + geom_line(data = df,
                               mapping = aes(x = Date, y = value, group = type,
                                   colour = type))
            g <- g + labs(colour = ylab)
            if (facets) {
                g <- g + facet_grid(type ~ ., scales = "free")
            }
        } else {
            g <- g + geom_line(data = df,
                               mapping = aes(x = Date, y = value, group = type),
                               alpha = alpha, colour = col1)
        }
        g <- g + labs(x = "", y = ylab)
    }
    g
   
}

## ************************************************************************************
##' @rdname autoplot.bts
##' @export
autolayer.bts <- function(object, col1 = "darkgray", ...) {
    
    value <- variable <- type <- Date <- NULL
    nc <- ncol(object)
    ylab <- attr(object, "label")
    if (is.null(ylab)) ylab <- ""
    col <- "orangered"
    alpha <- opacity(nc)
    
    if (nc == 1L) {
        df <- data.frame(Date = as.Date(attr(object, "date")),
                         value = unclass(object)[ , 1L])
        geom_line(data = df,
                  mapping = aes(x = Date, y = value),
                  alpha = alpha, colour = col)
    } else {
        
        if (FALSE) {
            df <- data.frame(date = as.Date(attr(object, "date")), object)
            df <- reshape2::melt(df, id.vars = "date")
        } else {
            df <- as.data.frame(object, longFormat = TRUE)
        }
        
        clab <- attr(object, "collabels")
        if (!is.null(clab)) levels(df$variable) <- clab
        
        if (nc <= 12L) {
            geom_line(data = df,
                      mapping = aes(x = Date, y = value, group = type,
                          colour = type))
        } else {
            geom_line(data = df,
                      mapping = aes(x = Date, y = value, group = type),
                      alpha = alpha, colour = col1)
        }
    }
   
}

## *****************************************************************************
##' Autoplot for Block Functional Time Series \code{bfts} objects.
##'
##' @title Autoplot for Block Functional Time Series \code{bfts} Objects
##' 
##' @param object A \code{bfts} object such as returned by the
##' \code{density} and \code{cdf} methods for \code{TVGEV} objects.
##'
##'
##' @param col Colour for the line and fill.
##'
##' @param fill If non \code{NA}, the region between the curve and the
##' axis will be filled as is quite usual with probability densities.
##' 
##' @param ... Not used yet.
##'
##' @return A \code{ggplot} object or nothing.
##' 
##' @note This function is mainly intended to be used for density
##' functions or cdf. In the former case, \code{fill = TRUE} may be
##' preferred.
##'
##' @export
##' 
autoplot.bfts <- function(object,
                          col = "orangered",
                          fill = NA, ...) {
    Fun <- Date <- x <- NULL
    
    if (nrow(object) > 12) {
        stop("too many rows")
    }
    d <- format(attr(object, "date"))
    ff <- data.frame(Date = rep(d, each = ncol(object)),
                     x = rep(attr(object, "x"), times = length(d)),
                     Fun = as.vector(t(object)))
    
    g <- ggplot(data = ff)
    g <- g + geom_line(mapping = aes(x = x, y = Fun, group = Date),
                       colour = "orangered")
    if (!is.na(fill)) {
        g <- g + geom_ribbon(mapping = aes(x = x, ymin = 0, ymax = Fun,
                                 group = Date),
                             fill = "orangered", alpha = 0.4)
    }
    g <- g + labs(y = attr(object, "label"), x = "")
    if (nrow(object) > 1L) g <- g + coord_flip()
    g <- g + facet_wrap( ~ Date)
    g
    
}

##' @method autolayer bfts
##' @export
autolayer.bfts <- function(object, col = "orangered", fill = NA, ...) {
    stop("The 'autolayer' method is not yet implemtented for the class ",
         "\"bfts\"")
}
