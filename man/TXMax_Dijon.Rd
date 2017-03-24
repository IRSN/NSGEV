\name{TXMax_Dijon}
\alias{TXMax_Dijon}
\docType{data}
\encoding{UTF8}
\title{
  
  Annual Maximum Temperature in Dijon 

}
\description{
  
  Annual Maxima for the series of daily Temperature maXimum (TX) in
  Dijon (Burgundy, France).

}
\usage{data("TXMax_Dijon")}
\format{

  A data frame with 96 observations on the following 2 variables.

  \describe{
    \item{\code{Year}}{

      Numeric vector giving the year.

    }
    \item{\code{TXMax}}{

      Numeric vector (with \code{NA}) giving the annual maximum of the
      daily maximum temperature (Celsius).

    }
  }
}
\details{

  The temperature in Celsius is recorded at Dijon-Longvic (France) one
  of the official french meteorological stations.

}
\source{
  
  The daily \code{TX} data are provided by European Data and Cimate
  Assessment \url{http://www.ecad.eu}, and were initially provided by
  Météo-France. They are available witin a \emph{Blended ECA dataset}
  for \emph{daily data}.

  The annual maxima were computed excluding the years with a missing
  daily observation. An exception is the year 2016 for which one
  observation was missing (2016-06-16): the missing TX was checked to be
  smaller than the maxima over non-missing observations
  
}
%% \references{}
\examples{
plot(TXMax ~ Year, data = TXMax_Dijon, type = "o", pch = 16,
     col = "orangered", main = "Annual Maxima of Temperature in Dijon")

## Use the begining of a year as timestamp.
df <- within(TXMax_Dijon, Date <- as.Date(sprintf("\%4d-01-01", Year)))
}
\keyword{datasets}
