\name{gw.adapt}
\alias{gw.adapt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Adaptive kernel for GWR}
\description{
  The function constructs weights using an adaptive kernal for geographically weighted regression
}
\usage{
gw.adapt(dp, fp, quant, longlat=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp}{data points coordinates}
  \item{fp}{fit points coordinates}
  \item{quant}{proportion of data points to include in the weights}
  \item{longlat}{if TRUE, use distances on an ellipse with WGS84 parameters}
}


\value{
  a vector of weights
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\keyword{spatial}

