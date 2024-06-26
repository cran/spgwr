## ----echo=FALSE, results="hide"-----------------------------------------------
load(system.file("backstore/nyGWR.RData", package="spgwr"))

## -----------------------------------------------------------------------------
library(spgwr)

## ----eval=FALSE---------------------------------------------------------------
#  if (packageVersion("spData") >= "2.3.2") {
#      NY8a <- sf::st_read(system.file("shapes/NY8_utm18.gpkg", package="spData"))
#  } else {
#      NY8a <- sf::st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData"))
#      sf::st_crs(NY8a) <- "EPSG:32618"
#      NY8a$Cases <- NY8a$TRACTCAS
#  }
#  NY8 <- as(NY8a, "Spatial")

## ----eval=FALSE---------------------------------------------------------------
#  bwG <- gwr.sel(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, gweight=gwr.Gauss, verbose=FALSE)
#  gwrG <- gwr(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, bandwidth=bwG, gweight=gwr.Gauss, hatmatrix=TRUE)

## -----------------------------------------------------------------------------
gwrG

## ----eval=FALSE---------------------------------------------------------------
#  gbwG <- ggwr.sel(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", gweight=gwr.Gauss, verbose=FALSE)
#  ggwrG <- ggwr(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", bandwidth=gbwG, gweight=gwr.Gauss)

## -----------------------------------------------------------------------------
ggwrG

## ----fig.caption="GWR local coefficient estimates for the exposure to TCE site covariate"----
spplot(ggwrG$SDF, "PEXPOSURE", col.regions=grey.colors(7, 0.95, 0.55, 2.2), cuts=6)

## ----fig.caption="plots of GWR local coefficient estimates showing the effects of GWR collinearity forcing"----
pairs(as(ggwrG$SDF, "data.frame")[,2:5])

