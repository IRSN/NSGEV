## *************************************************************************   
##' Print method for \code{bts} objects representing Block Time Series.
##'
##' While \code{bts} objects are inheriting from the \code{"matrix"}
##' class, their \code{print} method aims to show the date range in
##' compact format a dozen of lines with limited text width.
##'
##' @title Print Method for \code{bts} Objects Representing Block Time
##' Series
##' 
##' @param x A \code{bts} object.
##'
##' @param ... Not used yet.
##'
##' @return Nothing.
##' 
print.bts <- function(x, ...) {
    
    date <- format(attr(x, "date"))
    cat("Block time series\n\n")
    cat(sprintf("o First block: %s\no Last block:  %s\n\n",
                date[1], date[length(date)]))
    x <- noquote(format(x))
    
    skipCols <- FALSE
    if (ncol(x) > 12L) {
        skipCols <- TRUE
        indc1 <- 1:3
        indc2 <- ncol(x) - 2:0
    }
    
    skipRows <- FALSE
    indr1 <- 1L:nrow(x)
    if (nrow(x) > 6L) {
        skipRows <- TRUE
        indr1 <- 1:3L
        indr2 <- nrow(x) - 2:0
    } 
    
    ncc <- max(nchar(x))
    form <- paste("%", ncc, "s", collapse = "", sep = "")
    ## print(form)

    if (is.null(colnames(x))) colnames(x) <- letters[1:ncol(x)]
    if (skipCols) {
        cat(sprintf("%10s", "date"),
            sprintf(form, colnames(x)[indc1]),
            "   ...   ",
            sprintf(form, colnames(x)[indc2]), "\n")
        for (i in indr1) {
            cat(date[i],  x[i, indc1], "   ...   " , x[i, indc2], "\n")
        }
    } else {
        cat(sprintf("%10s", "date"),
            sprintf(form, colnames(x)), "\n")
        for (i in indr1) {
            cat(date[i], x[i, ], "\n")
        }
    }
    if (skipRows) {
        cat("  ... <more rows>\n")
        if (skipCols) {
            for (i in indr2) {
                cat(date[i],  x[i, indc1], "   ...   " , x[i, indc2], "\n")
            }
        } else {
            for (i in indr2) {
                cat(date[i], x[i, ], "\n")
            }
        }
    }
    
}


## *****************************************************************************
##' Time-plot of a \code{bts} object. 
##'
##' @title Time-Plot of \code{bts} Objects for Block Time Series
##'
##' @aliases plot.bts
##' @aliases lines.bts
##' 
##' @param x A \code{bts} object representing a Block Time Series,
##' such as the GEV parameters of a \code{TVGEV} model.
##' 
##' @param y Not used yet.
##'
##' @param gg Logical. If \code{TRUE}, a \code{ggplot} plot will be
##' drawn and returned by the function. If \code{FALSE} a standard
##' plot is shown and the function returns nothing.
##'
##' @param ... Further parameter passed to graphical functions.
##'
##' @return An object inheriting from \code{"ggplot"} if \code{gg}
##' is \code{TRUE}.
##'
##' @note When a plot is built with \code{gg = TRUE}, using the
##' \code{lines} method later will have no effect, because this method
##' is designed to work on standard graphics.
##'
##' @seealso \code{\link{quantile.TVGEV}}, \code{\link{mean.TVGEV}} for
##' some examples of \code{"bts"} objects and their use in plots.
##'
##' @examples
##' example(TVGEV)
##' plot(coef(res1, type = "theta"))
plot.bts <- function(x, y, gg = TRUE, col1 = "darkgray", ...) {

    value <- variable <- NULL
    nc <- ncol(x)
    ylab <- attr(x, "label")
    if (is.null(ylab)) ylab <- ""
    col <- "orangered"
    alpha <- opacity(nc)
    
    if (gg) {
        
        if (nc == 1L) {
            df <- data.frame(date = as.Date(attr(x, "date")), value = x[ , 1L])
            g <- ggplot(data = df, mapping = aes(x = date, y = value))
            g <- g + labs(colour = ylab, x = "", y = ylab)
            g <- g + geom_line(alpha = alpha, colour = col)
        } else {

            df <- data.frame(date = as.Date(attr(x, "date")), x)
            df <- melt(df, id.vars = "date")
            clab <- attr(x, "collabels")
            if (!is.null(clab)) levels(df$variable) <- clab
            
            if (nc <= 12L) {
                g <- ggplot(data = df,
                        mapping = aes(x = date, y = value, group = variable,
                            colour = variable))
                g <- g + labs(colour = ylab)
                g <- g + geom_line()
            } else {
                g <- ggplot(data = df,
                            mapping = aes(x = date, y = value, group = variable))
                g <- g + geom_line(alpha = alpha, colour = col1)
            }
            g <- g + labs(x = "", y = ylab)
        }
        g

    } else {
        date <- as.Date(attr(x, "date"))
        plot(range(date), range(x, na.rm = TRUE),
             type = "n",
             xlab = "", ylab = "quantile") 
        
        matlines(as.numeric(date), x, type = "l", lty = 1)
    }
        
}

## ****************************************************************************
##' @rdname plot.bts
##'
##' @param col1 Colour to be used when there are too many time series
##' to use one colour for each.
##' 
##' @param alpha Level of opacity. By default, the opacity is lowered
##' when the number of time series is increased.
##' 
lines.bts <- function(x, y, col1 = "gray", alpha = NULL, ...) {
    nc <- ncol(x)
    if (missing(alpha)) {
        alpha <- opacity(nc)
    }
    date <- as.Date(attr(x, "date"))
    if (nc <= 12) {
        matlines(as.numeric(date), x, type = "l", lty = 1,
                 ...)
    } else {
        col1 <- translude(col1, alpha = alpha)
        matlines(as.numeric(date), x, col = col1, type = "l", lty = 1,
                 ...)
    }
}

## *****************************************************************************
##' Plot for Block Functional Time Series \code{bfts} objects.
##'
##' @title Plot for Block Functional Time Series \code{bfts} Objects
##' 
##' @param x A \code{bfts} object such as returned by the \code{density} and
##' \code{cdf} methods for \code{TVGEV} objects.
##'
##' @param y Not used yet.
##'
##' @param gg Logical. If \code{TRUE}, a \code{ggplot} object is shown
##' and returned.
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
##' 
plot.bfts <- function(x, y, gg = TRUE,
                      col = "orangered",
                      fill = NA, ...) {
    Fun <- Date <- NULL
    
    if (!gg) {
        stop("For now only 'gg = TRUE' is possible.") 
    }
    
    if (nrow(x) > 12) {
        stop("too many rows")
    }
    d <- format(attr(x, "date"))
    ff <- data.frame(Date = rep(d, each = ncol(x)),
                     x = rep(attr(x, "x"), times = length(d)),
                     Fun = as.vector(t(x)))
    g <- ggplot(data = ff)
    g <- g + geom_line(mapping = aes(x = x, y = Fun, group = Date),
                       colour = "orangered")
    if (!is.na(fill)) {
        g <- g + geom_ribbon(mapping = aes(x = x, ymin = 0, ymax = Fun, group = Date),
                             fill = "orangered", alpha = 0.4)
    }
    g <- g + labs(y = attr(x, "label"), x = "")
    if (nrow(x) > 1L) g <- g + coord_flip()
    g <- g + facet_wrap( ~ Date)
    g
    
}
