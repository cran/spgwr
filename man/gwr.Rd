\name{gwr}
\alias{gwr}
\alias{print.gwr}
\title{Geographically weighted regression}
\description{
The function implements the basic geographically weighted regression approach
to exploring spatial non-stationarity for given global bandwidth and chosen
weighting scheme.
}
\usage{
gwr(formula, data=list(), coords, bandwidth, gweight=gwr.gauss, 
	adapt=NULL, hatmatrix = FALSE, fit.points, 
	longlat=FALSE, se.fit=FALSE, weights, cl=NULL)
\method{print}{gwr}(x, \dots)
}

\arguments{
  \item{formula}{regression model formula as in \code{lm}}
  \item{data}{model data frame, or SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations; may be omitted if the object passed through the data argument is from package \pkg{sp}}
  \item{bandwidth}{bandwidth used in the weighting function, possibly
  calculated by \code{gwr.sel}}
  \item{gweight}{geographical weighting function, at present only
  \code{gwr.gauss()} or \code{gwr.bisquare()}}
  \item{adapt}{either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours)}
  \item{hatmatrix}{if TRUE, return the hatmatrix as a component of the result}
  \item{fit.points}{an object containing the coordinates of fit points; often an object from package \pkg{sp}; if missing, the coordinates given through the data argument object, or the coords argument are used}
  \item{longlat}{if TRUE, use distances on an ellipse with WGS84 parameters}
  \item{se.fit}{if TRUE, return local coefficient standard errors}
  \item{weights}{case weights used as in weighted least squares, beware of scaling issues}
  \item{cl}{if NULL, ignored, otherwise \code{cl} must be an object describing a \dQuote{cluster} created using \code{makeCluster} in the snow package. The cluster will then be used to hand off the calculation of local coefficients to cluster nodes, if fit points have been given as an argument, and \code{hatmatrix=FALSE}}
  \item{x}{an object of class "gwr" returned by the \code{gwr} function}
  \item{...}{arguments to be passed to other functions}
}
\details{
The function applies the weighting function in turn to each of the 
observations, or fit points if given, calculating a weighted regression 
for each point. The results may be explored to see if coefficient values vary over space. The local coefficient estimates may be made on a multi-node cluster using the \code{cl} argument to pass through a \dQuote{snow} cluster. The function will then divide the fit points (which must be given separately) between the clusters for fitting. Note that each node will need to have the \dQuote{spgwr} package present, so initiating by \code{clusterEvalQ(cl, library(spgwr))} may save a little time per node. The function clears the global environment on the node of objects sent. Using two nodes reduces timings to a little over half the time for a single node.}

\value{
A list of class \dQuote{gwr}:
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package "sp") with fit.points, weights, GWR coefficient estimates, R-squared, and coefficient standard errors in its "data" slot.}
  \item{lhat}{Leung et al. L matrix}
  \item{lm}{Ordinary least squares global regression on the same model formula, as returned by lm.wfit().}
  \item{bandwidth}{the bandwidth used.}
  \item{this.call}{the function call used.}
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley;
\url{http://www.nuim.ie/ncg/GWR/index.htm}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr.sel}}, \code{\link{gwr.gauss}},
\code{\link{gwr.bisquare}}}

\examples{
data(columbus)
col.lm <- lm(crime ~ income + housing, data=columbus)
summary(col.lm)
col.bw <- gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y))
col.gauss <- gwr(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), bandwidth=col.bw, hatmatrix=TRUE)
col.gauss
col.d <- gwr.sel(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), gweight=gwr.bisquare)
col.bisq <- gwr(crime ~ income + housing, data=columbus,
  coords=cbind(columbus$x, columbus$y), bandwidth=col.d, 
  gweight=gwr.bisquare, hatmatrix=TRUE)
col.bisq
data(georgia)
g.adapt.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + PctBlack, data=gSRDF, adapt=TRUE)
res.adpt <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + PctBlack, data=gSRDF, adapt=g.adapt.gauss)
res.adpt
pairs(as(res.adpt$SDF, "data.frame")[,2:8], pch=".")
brks <- c(-0.25, 0, 0.01, 0.025, 0.075)
cols <- grey(5:2/6)
plot(res.adpt$SDF, col=cols[findInterval(res.adpt$SDF$PctBlack, brks, all.inside=TRUE)])
\dontrun{
  g.bw.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + PctBlack, data=gSRDF)
  res.bw <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + PctBlack, data=gSRDF, bandwidth=g.bw.gauss)
  res.bw
  pairs(as(res.bw$SDF, "data.frame")[,2:8], pch=".")
  plot(res.bw$SDF, col=cols[findInterval(res.bw$SDF$PctBlack, brks, all.inside=TRUE)])
  g.bw.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + PctBlack, data=gSRDF, longlat=TRUE)
  if (suppressWarnings(require(maptools)) &&
    suppressWarnings(require(gpclib))) {
    gSR <- as(gSRDF, "SpatialPolygons")
    length(getSpPpolygonsSlot(gSR))
    gSRouter <- unionSpatialPolygons(gSR, IDs=rep("Georgia", 159))
    SG <- GE_SpatialGrid(gSRouter, maxPixels = 100)
    SPxMASK0 <- overlay(gSRouter, SG$SG)
    SGDF <- SpatialGridDataFrame(slot(SG$SG, "grid"),
      data=data.frame(SPxMASK0=SPxMASK0),
      proj4string=CRS(proj4string(gSRouter)))
    SPxDF <- as(SGDF, "SpatialPixelsDataFrame")
    res.bw <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov +
      PctBlack, data=gSRDF, bandwidth=g.bw.gauss, fit.points=SPxDF,
      longlat=TRUE)
    res.bw
    spplot(res.bw$SDF, "PctBlack")
  }
}
}
\keyword{spatial}

