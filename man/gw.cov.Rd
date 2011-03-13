\name{gw.cov}
\alias{gw.cov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Geographically weighted local statistics}
\description{
  The function provides an implementation of geographically weighted local statistics based on Chapter 7 of the GWR book - see references. Local means, local standard deviations, local standard errors of the mean, standardised differences of the global and local means, and local covariances and if requested correlations, are reported for the chosed fixed or adaptive bandwidth and weighting function.
}
\usage{
gw.cov(x, vars, fp, adapt = NULL, bw, gweight = gwr.bisquare, cor = TRUE, var.term = FALSE, longlat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{x} should be a SpatialPolygonsDataFrame object or a SpatialPointsDataFrame object }
  \item{vars}{ \code{vars} is a vector of column numbers or a vector of column names applied to the columns of the data frame in the data slot of x }
  \item{fp}{ \code{fp} if given contains the fit points to be used, for example a SpatialPixels object describing the grid of points to be used }
  \item{adapt}{ \code{adapt} if given should lie between 0 and 1, and indicates the proportion of observations to be included in the weighted window - it cannot be selected automatically}
  \item{bw}{ \code{bw} when adapt is not given, the bandwidth chosen to suit the data set - it cannot be selected automatically}
  \item{gweight}{ \code{gweight} default gwr.bisquare - the weighting function to use }
  \item{cor}{ \code{cor} default TRUE, report correlations in addition to covariances }
  \item{var.term}{ \code{var.term} default FALSE, if TRUE apply a correction to the variance term }
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
}

\value{
  If argument fp is given, and it is a SpatialPixels object, a SpatialPixelsDataFrame is returned, if it is any other coordinate object, a SpatialPointsDataFrame is returned. If argument fp is not given, the object returned will be the class of object x. The data slot will contain a data frame with local means, local standard deviations, local standard errors of the mean, standardised differences of the global and local means, and local covariances and if requested correlations.
}

\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley (chapter 7);
\url{http://www.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr}}}

\examples{
data(georgia)
SRgwls <- gw.cov(gSRDF, vars=6:11, bw=2, longlat=FALSE)
names(SRgwls$SDF)
spplot(SRgwls$SDF, "mean.PctPov")
spplot(SRgwls$SDF, "sd.PctPov")
spplot(SRgwls$SDF, "sem.PctPov")
spplot(SRgwls$SDF, "diff.PctPov")
spplot(SRgwls$SDF, "cor.PctPov.PctBlack.")
SRgwls <- gw.cov(gSRDF, vars=6:11, bw=150, longlat=TRUE)
names(SRgwls$SDF)
spplot(SRgwls$SDF, "mean.PctPov")
spplot(SRgwls$SDF, "sd.PctPov")
spplot(SRgwls$SDF, "sem.PctPov")
spplot(SRgwls$SDF, "diff.PctPov")
spplot(SRgwls$SDF, "cor.PctPov.PctBlack.")
data(gSRouter)
gGrid <- sample.Polygons(slot(gSRouter, "polygons")[[1]], 5000,
  type="regular")
gridded(gGrid) <- TRUE
SGgwls <- gw.cov(gSRDF, vars=6:11, fp=gGrid, bw=150, longlat=TRUE)
names(SGgwls$SDF)
spplot(SGgwls$SDF, "mean.PctPov")
spplot(SGgwls$SDF, "sd.PctPov")
spplot(SGgwls$SDF, "sem.PctPov")
spplot(SGgwls$SDF, "diff.PctPov")
spplot(SGgwls$SDF, "cor.PctPov.PctBlack.")
}
\keyword{spatial}

