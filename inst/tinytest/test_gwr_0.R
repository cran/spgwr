library(spgwr)
data(columbus, package="spData")
expect_silent(col.bw <- gwr.sel(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), verbose = FALSE))
expect_true(isTRUE(all.equal(col.bw, 2.275059624)))
expect_silent(col.gauss <- gwr(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), bandwidth=col.bw, hatmatrix=TRUE))
expect_true(isTRUE(all.equal(col.gauss$results$edf, 19.38370134)))
expect_warning(col.d <- gwr.sel(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), gweight=gwr.bisquare, verbose = FALSE))
expect_true(isTRUE(all.equal(col.d, 33.07036513)))
expect_silent(col.bisq <- gwr(CRIME ~ INC + HOVAL, data=columbus,
  coords=cbind(columbus$X, columbus$Y), bandwidth=col.d, 
  gweight=gwr.bisquare, hatmatrix=TRUE))
expect_true(isTRUE(all.equal(col.bisq$results$edf, 44.39511877)))
