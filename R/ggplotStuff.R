##' Create a \code{ggplot} for a \code{quantMax.TVGEV} as created by
##' using the \code{quantMax} method.
##' 
##' @title Create a \code{ggplot} for a \code{quantMax.TVGEV} Object
##' 
##' @method autoplot quantMax.TVGEV
##' @export
##' 
##' @param object An object with class \code{"quantMax.TVGEV"} as
##'     produced by the \code{\link{quantMax}} method.
##' 
##' @param fillConf Logical. If \¢ode{TRUE} the "confidence band(s)"
##'     will be filled.
##'
##' @param ... Not used.
##'
##' @seealso The \code{\link{quantMax.TVGEV}} method.
##' 
autoplot.quantMax.TVGEV <- function(object, fillConf = TRUE, ...) {

    L <- U <- Quant <- ProbExc <- NULL
    
    ## Arrange the confidence levels so that the greatest level
    ## will be plotted first
    u <- sort(unique(object$Level), decreasing = TRUE)
    u <- formatLevel(u)
    if (length(u) > 2) {
        warning("a maximum of two confidence levels should be used")
    }
                
    object <- within(object,
                     Level <- factor(formatLevel(Level), levels = u))
    
    gg <- ggplot(data = object)

    if (!is.null(object$L) && !is.null(object$U)) {
        
        if (fillConf) {
            gg <- gg +
                geom_ribbon(mapping = aes(x = ProbExc, ymin = L, ymax = U,
                                          group = Level, fill = Level),
                            alpha = 0.2)
            ## gg <- gg + scale_fill_manual(values = c("SteelBlue1", "SteelBlue3"))
            gg <- gg + scale_fill_manual(name = "Level",
                                         values = c("gray75", "gray60"))
            
        }

        gg <- gg +
            geom_line(mapping = aes(x = ProbExc, y = L, linetype = Level,
                                       group = Level, colour = Level),
                      alpha = 0.9, size = 0.8)
        gg <- gg +
            geom_line(mapping = aes(x = ProbExc, y = U, linetype = Level,
                                    group = Level, colour = Level),
                      alpha = 0.9, linewidth = 0.8)
        
        gg <- gg + scale_colour_manual(values = c("gray75", "gray60"))
        
    }
    
    gg <- gg + scale_x_continuous(trans = .gumbel_trans_p,
                                  breaks = .gumBreaks_p,
                                  minor_breaks = .gumBreaks_p)
    gg <- gg + geom_line(mapping = aes(x = ProbExc, y = Quant),
                         col = "orangered")

    gg <- gg + xlab("Prob. of exceedance") + ylab("Quantile") 
    
}
