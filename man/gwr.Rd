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
gwr(formula, data=list(), coords, bandwidth, gweight=gwr.Gauss, 
	adapt=NULL, hatmatrix = FALSE, fit.points, longlat=NULL, 
	se.fit=FALSE, weights, cl=NULL, predictions = FALSE, 
        fittedGWRobject = NULL, se.fit.CCT = TRUE)
\method{print}{gwr}(x, \dots)
}

\arguments{
  \item{formula}{regression model formula as in \code{lm}}
  \item{data}{model data frame, or SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package \pkg{sp}}
  \item{coords}{matrix of coordinates of points representing the spatial
  positions of the observations; may be omitted if the object passed through the data argument is from package \pkg{sp}}
  \item{bandwidth}{bandwidth used in the weighting function, possibly
  calculated by \code{gwr.sel}}
  \item{gweight}{geographical weighting function, at present 
  \code{gwr.Gauss()} default, or \code{gwr.gauss()}, the previous default or \code{gwr.bisquare()}}
  \item{adapt}{either NULL (default) or a proportion between 0 and 1 of observations to include in weighting scheme (k-nearest neighbours)}
  \item{hatmatrix}{if TRUE, return the hatmatrix as a component of the result, ignored if \code{fit.points} given}
  \item{fit.points}{an object containing the coordinates of fit points; often an object from package \pkg{sp}; if missing, the coordinates given through the data argument object, or the coords argument are used}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers; if x is a SpatialPoints object, the value is taken from the object itself}
  \item{se.fit}{if TRUE, return local coefficient standard errors - if \code{hatmatrix} is TRUE and no fit.points are given, two effective degrees of freedom sigmas will be used to generate alternative coefficient standard errors}
  \item{weights}{case weights used as in weighted least squares, beware of scaling issues, probably unsafe}
  \item{cl}{if NULL, ignored, otherwise \code{cl} must be an object describing a \dQuote{cluster} created using \code{makeCluster} in the \pkg{parallel} package. The cluster will then be used to hand off the calculation of local coefficients to cluster nodes, if fit points have been given as an argument, and \code{hatmatrix=FALSE}}
  \item{predictions}{default FALSE; if TRUE and no fit points given, return GW fitted values at data points, if fit points given and are a Spatial*DataFrame object containing the RHS variables in the formula, return GW predictions at the fit points}
  \item{fittedGWRobject}{a fitted \code{gwr} object with a hatmatrix (optional), if given, and if fit.points are given and if \code{se.fit} is TRUE, two effective degrees of freedom sigmas will be used to generate alternative coefficient standard errors}
  \item{se.fit.CCT}{default TRUE, compute local coefficient standard errors using formula (2.14), p. 55, in the GWR book}
  \item{x}{an object of class "gwr" returned by the \code{gwr} function}
  \item{...}{arguments to be passed to other functions}
}
\details{
The function applies the weighting function in turn to each of the 
observations, or fit points if given, calculating a weighted regression 
for each point. The results may be explored to see if coefficient values vary over space. The local coefficient estimates may be made on a multi-node cluster using the \code{cl} argument to pass through a \pkg{parallel} cluster. The function will then divide the fit points (which must be given separately) between the clusters for fitting. Note that each node will need to have the \dQuote{spgwr} package present, so initiating by \code{clusterEvalQ(cl, library(spgwr))} may save a little time per node. The function clears the global environment on the node of objects sent. Using two nodes reduces timings to a little over half the time for a single node.

The section of the examples code now includes two simulation scenarios, showing how important it is to check that mapped pattern in local coefficients is actually there, rather than being an artefact.
}

\value{
A list of class \dQuote{gwr}:
  \item{SDF}{a SpatialPointsDataFrame (may be gridded) or SpatialPolygonsDataFrame object (see package "sp") with fit.points, weights, GWR coefficient estimates, R-squared, and coefficient standard errors in its "data" slot.}
  \item{lhat}{Leung et al. L matrix}
  \item{lm}{Ordinary least squares global regression on the same model formula, as returned by lm.wfit().}
  \item{bandwidth}{the bandwidth used.}
  \item{this.call}{the function call used.}
}
\references{Fotheringham, A.S., Brunsdon, C., and Charlton, M.E., 2002,
Geographically Weighted Regression, Chichester: Wiley; Paez A, Farber S, Wheeler D, 2011, "A simulation-based study of geographically weighted regression as a method for investigating spatially varying relationships", Environment and
Planning A 43(12) 2992-3010;
\url{http://gwr.nuim.ie/}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{\code{\link{gwr.sel}}, \code{\link{gwr.gauss}},
\code{\link{gwr.bisquare}}}

\examples{
data(columbus, package="spData")
col.lm <- lm(CRIME ~ INC + HOVAL, data=columbus)
summary(col.lm)
col.bw <- gwr.sel(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y))
col.gauss <- gwr(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), bandwidth=col.bw, hatmatrix=TRUE)
col.gauss
col.d <- gwr.sel(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), gweight=gwr.bisquare)
col.bisq <- gwr(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), bandwidth=col.d, 
  gweight=gwr.bisquare, hatmatrix=TRUE)
col.bisq
data(georgia)
g.adapt.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + 
 PctPov + PctBlack, data=gSRDF, adapt=TRUE)
res.adpt <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov + 
 PctBlack, data=gSRDF, adapt=g.adapt.gauss)
res.adpt
pairs(as(res.adpt$SDF, "data.frame")[,2:8], pch=".")
brks <- c(-0.25, 0, 0.01, 0.025, 0.075)
cols <- grey(5:2/6)
plot(res.adpt$SDF, col=cols[findInterval(res.adpt$SDF$PctBlack, brks,
 all.inside=TRUE)])
\donttest{
# simulation scenario with patterned dependent variable
set.seed(1)
X0 <- runif(nrow(gSRDF)*3)
X1 <- matrix(sample(X0), ncol=3)
X1 <- prcomp(X1, center=FALSE, scale.=FALSE)$x
gSRDF$X1 <- X1[,1]
gSRDF$X2 <- X1[,2]
gSRDF$X3 <- X1[,3]
bw <- gwr.sel(PctBach ~ X1 + X2 + X3, data=gSRDF, verbose=FALSE)
out <- gwr(PctBach ~ X1 + X2 + X3, data=gSRDF, bandwidth=bw, hatmatrix=TRUE)
out
spplot(gSRDF, "PctBach", col.regions=grey.colors(20))
spplot(gSRDF, c("X1", "X2", "X3"), col.regions=grey.colors(20))
# pattern in the local coefficients
spplot(out$SDF, c("X1", "X2", "X3"), col.regions=grey.colors(20))
# but no "significant" pattern
spplot(out$SDF, c("X1_se", "X2_se", "X3_se"), col.regions=grey.colors(20))
out$SDF$X1_t <- out$SDF$X1/out$SDF$X1_se
out$SDF$X2_t <- out$SDF$X2/out$SDF$X2_se
out$SDF$X3_t <- out$SDF$X3/out$SDF$X3_se
spplot(out$SDF, c("X1_t", "X2_t", "X3_t"), col.regions=grey.colors(20))
# simulation scenario with random dependent variable
yrn <- rnorm(nrow(gSRDF))
gSRDF$yrn <- sample(yrn)
bw <- gwr.sel(yrn ~ X1 + X2 + X3, data=gSRDF, verbose=FALSE)
# bandwidth selection maxes out at 620 km, equal to upper bound
# of line search
out <- gwr(yrn ~ X1 + X2 + X3, data=gSRDF, bandwidth=bw, hatmatrix=TRUE)
out
spplot(gSRDF, "yrn", col.regions=grey.colors(20))
spplot(gSRDF, c("X1", "X2", "X3"), col.regions=grey.colors(20))
# pattern in the local coefficients
spplot(out$SDF, c("X1", "X2", "X3"), col.regions=grey.colors(20))
# but no "significant" pattern
spplot(out$SDF, c("X1_se", "X2_se", "X3_se"), col.regions=grey.colors(20))
out$SDF$X1_t <- out$SDF$X1/out$SDF$X1_se
out$SDF$X2_t <- out$SDF$X2/out$SDF$X2_se
out$SDF$X3_t <- out$SDF$X3/out$SDF$X3_se
spplot(out$SDF, c("X1_t", "X2_t", "X3_t"), col.regions=grey.colors(20))
# end of simulations
}
\donttest{
data(meuse)
coordinates(meuse) <- c("x", "y")
meuse$ffreq <- factor(meuse$ffreq)
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
meuse.grid$ffreq <- factor(meuse.grid$ffreq)
gridded(meuse.grid) <- TRUE
xx <- gwr(cadmium ~ dist, meuse, bandwidth = 228, hatmatrix=TRUE)
xx
x <- gwr(cadmium ~ dist, meuse, bandwidth = 228, fit.points = meuse.grid,
 predict=TRUE, se.fit=TRUE, fittedGWRobject=xx)
x
spplot(x$SDF, "pred")
spplot(x$SDF, "pred.se")
}
\dontrun{
  g.bw.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB +
    PctPov + PctBlack, data=gSRDF)
  res.bw <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov +
    PctBlack, data=gSRDF, bandwidth=g.bw.gauss)
  res.bw
  pairs(as(res.bw$SDF, "data.frame")[,2:8], pch=".")
  plot(res.bw$SDF, col=cols[findInterval(res.bw$SDF$PctBlack, brks,
    all.inside=TRUE)])
  g.bw.gauss <- gwr.sel(PctBach ~ TotPop90 + PctRural + PctEld + PctFB +
    PctPov + PctBlack, data=gSRDF, longlat=TRUE)
  data(gSRouter)
  require(maptools)
  SG <- GE_SpatialGrid(gSRouter, maxPixels = 100)
  SPxMASK0 <- over(SG$SG, gSRouter)
  SGDF <- SpatialGridDataFrame(slot(SG$SG, "grid"),
    data=data.frame(SPxMASK0=SPxMASK0),
    proj4string=CRS(proj4string(gSRouter)))
  SPxDF <- as(SGDF, "SpatialPixelsDataFrame")
  res.bw <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov +
    PctBlack, data=gSRDF, bandwidth=g.bw.gauss, fit.points=SPxDF,
    longlat=TRUE)
  res.bw
  res.bw$timings
  spplot(res.bw$SDF, "PctBlack")
  require(parallel)
  cl <- makeCluster(detectCores())
  res.bwc <- gwr(PctBach ~ TotPop90 + PctRural + PctEld + PctFB + PctPov +
    PctBlack, data=gSRDF, bandwidth=g.bw.gauss, fit.points=SPxDF,
    longlat=TRUE, cl=cl)
  res.bwc
  res.bwc$timings
  stopCluster(cl)
}
}
\keyword{spatial}

