

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
