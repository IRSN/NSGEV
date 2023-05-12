.gumBreaks_p <- c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01,  0.005, 0.002, 0.001)

##' @importFrom scales trans_new
.gumbel_trans_p <- scales::trans_new(name = "gumbel",
                                     transform = function(p) -log(-log(1 - p)),
                                     inverse = function(y) 1 - exp(-exp(-y)),
                                     breaks = myBreaks <- c(0.5, 0.2, 0.1,
                                                            0.05, 0.02, 0.01,  0.005, 0.002, 0.001),
                                     domain = c(0, 1))

.gumBreaks_m <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
.gumbel_trans_m <- scales::trans_new(name = "gumbel",
                                     transform = function(m) -log(-log(1 - 1 / m)),
                                     inverse = function(y) 1 / (1 - exp(-exp(-y))),
                                     breaks = myBreaks <- c(2, 5, 10,
                                                            20, 50, 100, 200, 500, 1000),
                                     domain = c(1, 1000))
