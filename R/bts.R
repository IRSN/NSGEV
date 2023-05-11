## **************************************************************************
##' Create a \code{bts} object representing a block time-series.
##'
##' This structure is intended to represent \emph{slowly varying}
##' time-series only. So we can plot the column of \code{bts} object
##' against the date vector given as attribute. Since a date
##' represents the beginning of a block, each polyline (or "curve")
##' actually shifted by an half-block.
##' 
##' @title Create a \code{bts} Object Representing a Block Time-Series
##' 
##' @param x A matrix with the time-series components as its columns.
##'
##' @param dateFrom A beginning date used when \code{date} is
##' missing. Then a date vector 
##'
##' @param date A vector that can be coerced to the class
##' \code{"Date"} with length \code{nrow(x)}.
##' 
##' @return An object with class \code{"bts"} inheriting from
##' \code{"matrix"}. This object has a \code{"date"} attribute which
##' is used to identify the beginning of the blocks.
##'
##' @note Some matrix operations do not make sense for \code{bts}
##' objects. A \code{"bts"} is transformed in an ordinary numeric
##' matrix by simply using \code{unclass}, see examples.
##'
##' @export
##' 
##' @examples
##' x <- matrix(1:300, ncol = 10)
##' myBts1 <- bts(x, dateFrom = "2010-01-01")
##' ## plot method
##' plot(myBts1)
##' mat <- unclass(myBts1)
##' class(mat)
##' x <- matrix(150 + 150 * runif(300), ncol = 100)
##' ## with more than 12 columns, the legend is no longer shown
##' myBts2 <- bts(x, date = c("2040-01-01", "2041-01-01", "2042-01-01"))
##' plot(myBts2)
bts <- function(x, dateFrom, date) {

    x <- as.matrix(x)
    nd <- nrow(x)

    if (missing(dateFrom)) {
        if (missing(date)) {

            if (!is.null(rownames(x))) {
                chck <- try(date <- as.Date(rownames(x)))
            }
            if (inherits(chck, "try-error")) {   
                stop("'x' does not have date rownames, so one of the two ",
                     "arguments 'date' and 'dateFrom' must be given")
            }
        }
        if (length(date) != nd) {
            stop("The length of 'date' must match the number of rows of 'x'") 
        }
        date <- as.Date(date)
    } else {
        if (!missing(date)) {
            warning("Since 'date' is given, 'dateFrom' is ignored")
        }
        date <- seq(from = as.Date(dateFrom), length.out = nd, by = "years")
    }
    
    attr(x, "date") <- date
    class(x) <- c("bts", "matrix")
    x

}


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
##' @method print bts
##' @export
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
##' Coerce an \code{bts} object into a data frame
##'
##' @title Coerce an \code{bts} Object into a Data Frame
##' 
##' @param x An object with class \code{bts} representing a by-block time
##' series.
##'
##' @param row.names Not used.
##' 
##' @param optional Not used.
##'
##' @param longFormat Logical. With the default \code{FALSE}, the
##' result has the same columns as \code{x} plus the new column
##' \code{Date}. When \code{longFormat} is \code{TRUE}, the data frame
##' is in "long format" with all the columns of \code{x} given in same
##' the column \code{value} of the result. A column \code{type} indicates
##' the name of the original column in \code{x}.
##' 
##' @param ... Not used yet.
##' 
##' @return A dataframe with a \code{Date} column containing the value
##' of the \code{"date"} attribute of \code{x}.
##'
##' @method as.data.frame bts
##' @export
##' 
##' @examples
##' myDate <- seq(from = as.Date("1996-01-01"),
##'             to = as.Date("2016-12-31"), by = "years")
##' X <- polynomX(myDate)
##' df <- as.data.frame(X)
##' head(df)
as.data.frame.bts <- function(x, row.names = NULL,
                              optional = FALSE,
                              longFormat = FALSE,
                                   ...) {

    type <- value <- Date <- NULL
    if ("Date" %in% colnames(x)) {
        stop("Invalid colnames \"Date\" in 'x'")
    }
    
    df <- as.data.frame.matrix(x)
    df <- data.frame(Date = as.Date(attr(x, "date")), df)
    rownames(df) <- NULL

    if (longFormat) {
        df <- tidyr::gather(df, key = type, value = value, -Date)
    }
    
    df
    
}

## ***************************************************************************************
##'  Coercion into a \code{bts} object.
##' 
##' @title Coercion into a \code{bts} Object
##'
##' @param object An object to be coerced into a \code{bts} object.
##'
##' @param ... Extra arguments for methods.
##'
##' @return An object with class \code{"bts"}.
##'
##' @export
##' 
as.bts <- function(object, ...) {
    UseMethod("as.bts")
}

## *************************************************************************************
##' Coerce a data frame into a \code{bts} object, using a suitable
##' column as a date.
##'
##' @title Coerce a Data Frame into a \code{bts} object
##'
##' @param object A data frame with one column
##'
##' @param dateName Character. Optional column name for a date in an
##' understandable format.
##'
##' @param yearName Character. Optional column name for a year. When
##' given, the column must be either numeric or character.
##'
##' @param ... Not used.
##'
##' @return An object with class \code{"bts"} inheriting from
##' \code{"matrix"}.
##'
##' @method as.bts data.frame
##' @export
##' 
##' @examples
##' ## does not work
##' bts0 <- try(as.bts(TXMax_Dijon))
##' ## works
##' bts1 <- try(as.bts(TXMax_Dijon, year = "Year"))
##' 
as.bts.data.frame <- function(object, dateName, yearName, ...) {
    
    if (!missing(dateName)) {
        
        if (!missing(yearName)) {
            stop("only one of the arguments 'dateName' and 'yearNames' ",
                 "can be given") 
        }
    
        ind <- match(dateName, colnames(object))
        if (length(ind) != 1) {
            stop("'dateName = ' must specify one column") 
        }
        .date <- try(as.Date(object[ , ind]))
        if (inherits(.date, "try-error")) {
            stop("Bad dateName")
        }
       
    } else if (!missing(yearName)) {
        ind <- match(yearName, colnames(object))
        if (length(ind) != 1) {
            stop("'yearName = ' must specify one column") 
        }
        if (is.character(object[ , ind])) {
            .date <- as.Date(sprintf("%s-01-01", object[ , ind]))
        } else if (is.numeric(object[ , ind])) {
            .date <- as.Date(sprintf("%4d-01-01", object[ , ind]))
        } else {
            stop("'year' can only specify a character of a ",
                 "numeric column")
        }
    } else {
        ind <- sapply(object, function(x) class(x)[1] == "Date")
        ind <- (1:ncol(object))[ind]
        if (length(ind) != 1) {
            stop("'dateName' an 'yearName' can both be missing only if 'object' ",
                 "contains exactly one column with class \"Date\"") 
        }
        .date <- object[ , ind]
    }
  
    x <- as.matrix(object[ , -ind, drop = FALSE])
    attr(x, "date") <- .date
    class(x) <- c("bts", "matrix")
    x

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
##' @param col1 A colour which is used in the case where \code{x} has
##' more than 12 columns. In this case, all the curves are drawn using
##' the colour \code{col1}.
##' 
##' @param facets Logical. If \code{True}, the columns of \code{x} will
##' be shown each in a facet. Only works if the number of columns is
##' between 2 and 12.
##'
##' @param ... Further parameter passed to graphical functions.
##'
##' @return An object inheriting from \code{"ggplot"} if \code{gg}
##' is \code{TRUE}.
##'
##' @note When a plot is built with \code{gg = TRUE}, using the
##' \code{lines} method later will have no effect, because this method
##' is designed to work on standard graphics, and not with graphics
##' produced by \pkg{ggplot2}.
##'
##' @seealso \code{\link{quantile.TVGEV}}, \code{\link{mean.TVGEV}} for
##' some examples of \code{"bts"} objects and their use in plots.
##'
##' @method plot bts
##' @export
##' 
##' @examples
##' example(TVGEV)
##' plot(coef(res1, type = "theta"))
plot.bts <- function(x, y, gg = TRUE, col1 = "darkgray", facets = FALSE, ...) {
    
    value <- variable <- Date <- NULL
    nc <- ncol(x)
    ylab <- attr(x, "label")
    if (is.null(ylab)) ylab <- ""
    col <- "orangered"
    alpha <- opacity(nc)
    
    if (gg) {

        if (((nc == 1L) || (nc > 12L)) && facets) {
            warning("Value of 'facets' ignored because the number of ",
                    "columns is 1 or is > 12")
        }
        
        if (nc == 1L) {
            df <- data.frame(Date = as.Date(attr(x, "date")), value = x[ , 1L])
            g <- ggplot(data = df, mapping = aes(x = Date, y = value))
            g <- g + labs(colour = ylab, x = "", y = ylab)
            g <- g + geom_line(alpha = alpha, colour = col)
        } else {

            if (FALSE) {
                df <- data.frame(date = as.Date(attr(x, "date")), x)
                df <- melt(df, id.vars = "date")
            } else {
                df <- as.data.frame(x)
                df <- reshape2::melt(df, id.vars = "Date")
            }
            
            clab <- attr(x, "collabels")
            if (!is.null(clab)) levels(df$variable) <- clab
            
            if (nc <= 12L) {
                g <- ggplot(data = df,
                        mapping = aes(x = Date, y = value, group = variable,
                            colour = variable))
                g <- g + labs(colour = ylab)
                g <- g + geom_line()
                if (facets) {
                    g <- g + facet_grid(variable ~ ., scales = "free")
                }
            } else {
                g <- ggplot(data = df,
                            mapping = aes(x = Date, y = value, group = variable))
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
##' @method plot bfts
##' @export
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
        g <- g + geom_ribbon(mapping = aes(x = x, ymin = 0, ymax = Fun,
                                 group = Date),
                             fill = "orangered", alpha = 0.4)
    }
    g <- g + labs(y = attr(x, "label"), x = "")
    if (nrow(x) > 1L) g <- g + coord_flip()
    g <- g + facet_wrap( ~ Date)
    g
    
}

##' Extracts parts of a \code{bts} object.
##'
##' @title Extract Parts of a \code{bts} Object
##'
##' @param x A \code{bts} object.
##' 
##' @param ... Elements to extract or replace. 
##'
##' @param drop \code{Logical}. If \code{TRUE} the result is coerced to the
##' lowest possible dimension.
##'
##' @method [ bts
##' @export
##' 
`[.bts` <- function(x, ..., drop = FALSE) {
    L <- list()
    attrs <- list()
    for (nm in c("date", "collabels", "label")) {
        attrs[[nm]] <- attr(x , which = nm)
    }
    x <- unclass(x)[..., drop = drop]
    ## x <- bts(x, date = attrs[["date"]])
    for (nm in c("date", "collabels", "label")) {
        attr(x, nm)  <- attrs[[nm]]
    }
    class(x) <- c("bts", "matrix")
    x
}
## ****************************************************************************
##' Extracts a subset of a \code{bts} object \code{x} observed between the
##' dates \code{start} and \code{end}.
##'
##' @title Extracts a Subset of a \code{bts} Object Observed
##' Between the Dates \code{start} and \code{end}
##'
##' @param x A \code{bts} object.
##'
##' @param start,end Objects coerced into a \code{Date} giving the start
##' and end.
##'
##' @param extend Logical.  If \code{TRUE}, the ‘start’ and ‘end’
##' values are allowed to extend the series.  If \code{TRUE}, attempts
##' to extend the series give an error.
##'
##' @param ... Not used
##'
##' @return A \code{bts} object.
##'
##' @importFrom stats window
##' @method window bts
##' @export
##' 
window.bts <- function(x, start, end, extend = FALSE, ... ) {

    date <- as.Date(attr(x, "date"))
    nd <- length(date)
    
    if (!missing(start)) {
        start <- as.Date(start)
    } else {
        start <- date[1]
    }

    if (!missing(end)) {
        end <- as.Date(end)
    } else {
        end <-  date[nd]
    }
        
    if (start > end) stop ("'start' must be >= 'end'")
    
    dateNew <- seq(from = start, to = end, by = "years")
    
    ind <- (date >= start) & (date <= end)
    xNew <- unclass(x)[ind, , drop = FALSE]
    if (start < date[1]) {
        if (extend) {
            xNew <- rbind(matrix(NA,
                                 nrow = length(dateNew[dateNew < date[1]]),
                                 ncol = ncol(xNew)), xNew)
        } else {
            stop("When 'extend' is FALSE, 'start' can not be before the ",
                 "start of the bts object")
        }
    } else if (start > date[1]) {
        dateNew <- dateNew[dateNew >= start]
    }
        
    if (end > date[nd]) {
        
        if (extend) {
            xNew <- rbind(xNew, matrix(NA,
                                       nrow = length(dateNew[dateNew > date[nd]]),
                                       ncol = ncol(xNew)))
        } else {
            stop("When 'extend' is FALSE, 'end' can not be after the ",
                 "end of the bts object")
        }
    } else if (end < date[nd]) {
        dateNew <- dateNew[dateNew <= end]
    }
        
    for (nm in c("collabels", "label")) {
        attr(xNew, nm)  <- attr(x , which = nm)
    }
    
    attr(xNew, "date") <- rownames(xNew) <- format(dateNew, "%Y-%m-%d")
    class(xNew) <- c("bts", "matrix")
    xNew

}

## ****************************************************************************
##' Coerce a \code{bts} object into a \code{ts} time series object.
##'
##' @title Coerce a \code{bts} Object into a \code{ts} Time Series Object
##'
##' @param x A \code{bts} object.
##'
##' @param ... Not used.
##' 
##' @return An object inheriting from the \code{"ts"} S3 class.
##'
##' @importFrom stats as.ts
##' @method as.ts bts
##' @export
##' 
as.ts.bts <- function(x, ...) {
    nc <- ncol(x)
    d <- as.Date(attr(x, "date"))
    startYear <- as.numeric(format(d[1], "%Y"))
    bd <- blockDuration(d)
    ts(data = x, start = startYear, frequency = 1 / bd)
}
