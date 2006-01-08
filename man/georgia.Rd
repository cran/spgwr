\name{georgia}
\alias{georgia}
\alias{gSRDF}
\docType{data}
\title{Georgia census data set (SpatialDataFramePolygons)}
\description{
  The Georgia census data set from Fotheringham et al. (2002) in shapefile format.
}
\usage{data(georgia)}
\format{
 A SpatialPolygonsDataFrame object (proj4string set to "+proj=latlong +datum=NAD27").

  The "data" slot is a data frame with 159 observations on the following 13 variables.
  \describe{
    \item{AreaKey}{a numeric vector}
    \item{Latitude}{a numeric vector}
    \item{Longitud}{a numeric vector}
    \item{TotPop90}{a numeric vector}
    \item{PctRural}{a numeric vector}
    \item{PctBach}{a numeric vector}
    \item{PctEld}{a numeric vector}
    \item{PctFB}{a numeric vector}
    \item{PctPov}{a numeric vector}
    \item{PctBlack}{a numeric vector}
    \item{ID}{a numeric vector}
    \item{X}{a numeric vector}
    \item{Y}{a numeric vector}
  }
}
\details{
Variables are from GWR3 file GeorgiaData.csv.
}
\source{
  \url{http://www.census.gov/geo/cob/bdy/co/co90shp/co13_d90_shp.zip}; \url{http://www.nuim.ie/ncg/GWR/index.htm}
}
\references{
  Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002, Geographically Weighted Regression: The Analysis of Spatially Varying Relationships, Chichester: Wiley.
}
\examples{
data(georgia)
plot(gSRDF)
}
\keyword{datasets}
