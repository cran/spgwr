# Copyright 2001-2007 Roger Bivand and Danlin Yu
# 

gwr.sel <- function(formula, data = list(), coords, adapt=FALSE, 
	gweight=gwr.gauss, method="cv", verbose=TRUE, longlat=FALSE) {
	if (!is.logical(adapt)) stop("adapt must be logical")
	if (!is.logical(longlat)) stop("longlat must be logical")
	if (is(data, "Spatial")) {
		if (!missing(coords))
		    warning("data is Spatial* object, ignoring coords argument")
		coords <- coordinates(data)
		data <- as(data, "data.frame")
	}
	if (missing(coords))
		stop("Observation coordinates have to be given")
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, method="model.frame", na.action=na.fail)
#	dist2 <- (as.matrix(dist(coords)))^2
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)
#	if (NROW(x) != NROW(dist2))
#		stop("Input data and coordinates have different dimensions")
	if (!adapt) {
		bbox <- cbind(range(coords[,1]), range(coords[,2]))
		difmin <- spDistsN1(bbox, bbox[2,], longlat)[1]
		if (any(!is.finite(difmin)))
			difmin[which(!is.finite(difmin))] <- 0
		beta1 <- difmin/1000
		beta2 <- difmin
		if (method == "cv") {
			opt <- optimize(gwr.cv.f, lower=beta1, upper=beta2, 
				maximum=FALSE, y=y, x=x, coords=coords, 
				gweight=gweight, verbose=verbose, longlat=longlat)
		} else {
			opt <- optimize(gwr.aic.f, lower=beta1, upper=beta2, 
				maximum=FALSE, y=y, x=x, coords=coords, 
				gweight=gweight, verbose=verbose, longlat=longlat)
		}
		bdwt <- opt$minimum
		res <- bdwt
	} else {
		beta1 <- 0
		beta2 <- 1
		if (method == "cv") {
			opt <- optimize(gwr.cv.adapt.f, lower=beta1, 
				upper=beta2, maximum=FALSE, y=y, x=x, 
				coords=coords, gweight=gweight, 
				verbose=verbose, longlat=longlat)
		} else {
			opt <- optimize(gwr.aic.adapt.f, lower=beta1, 
				upper=beta2, maximum=FALSE, y=y, x=x, 
				coords=coords, gweight=gweight, 
				verbose=verbose, longlat=longlat)
		}
		q <- opt$minimum
		res <- q
	}
	res
}

gwr.aic.f <- function(bandwidth, y, x, coords, gweight, verbose=TRUE, longlat=FALSE) {
    n <- NROW(x)
#    m <- NCOL(x)
    lhat <- matrix(nrow=n, ncol=n)
    flag <- 0
    options(show.error.messages = FALSE)
    for (i in 1:n) {
#        xx <- x[i, ]
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- 0
	w.i <- gweight(dxs^2, bandwidth)
#	w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
        lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
        if(!inherits(lm.i, "try-error")) {
            p <- lm.i$rank
	    p1 <- 1:p
	    inv.Z <- chol2inv(lm.i$qr$qr[p1, p1, drop=FALSE])
	    lhat[i,] <- t(x[i,]) %*% inv.Z %*% t(x) %*% diag(w.i)
        } else {
	    flag <- 1
	}
    }
    if (flag == 0) {
    	v1 <- sum(diag(lhat))
    	B1 <- t(diag(n)-lhat)%*%(diag(n)-lhat)
    	rss <- c(t(y)%*%B1%*%y)
    	sigma2.b <- rss / n
# NOTE 2* and sqrt() inserted for legibility
    	score <- 2*n*log(sqrt(sigma2.b)) + n*log(2*pi) + 
	    (n * (n + v1) / (n - 2 - v1))
    } else {
	score <- as.numeric(NA)
    }
    options(show.error.messages = TRUE)
    if (verbose) cat("Bandwidth:", bandwidth, "AIC:", score, "\n")
    score
}

gwr.cv.f <- function(bandwidth, y, x, coords, gweight, verbose=TRUE, longlat=FALSE)
{
    n <- NROW(x)
#    m <- NCOL(x)
    cv <- numeric(n)
    options(show.error.messages = FALSE)
    for (i in 1:n) {
        xx <- x[i, ]
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- 0
	w.i <- gweight(dxs^2, bandwidth)
#	w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bandwidth)
        w.i[i] <- 0
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
        lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
        if(!inherits(lm.i, "try-error")) {
            b <- coefficients(lm.i)
            cv[i] <- y[i] - (t(b) %*% xx)
        }
    }
    score <- sqrt(sum(t(cv) %*% cv)/n)
    options(show.error.messages = TRUE)
    if (verbose) cat("Bandwidth:", bandwidth, "CV score:", score, "\n")
    score
}

gwr.aic.adapt.f <- function(q, y, x, coords, gweight, verbose=TRUE, longlat=FALSE) {
    n <- NROW(x)
#    m <- NCOL(x)
    lhat <- matrix(nrow=n, ncol=n)
    bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
    flag <- 0
    options(show.error.messages = FALSE)
    for (i in 1:n) {
#        xx <- x[i, ]
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- 0
	w.i <- gweight(dxs^2, bw[i])
#	w.i <- gweight(spDistsN1(coords, coords[i,], longlat=longlat)^2, bw[i])
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
        lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
        if(!inherits(lm.i, "try-error")) {
            p <- lm.i$rank
	    p1 <- 1:p
	    inv.Z <- chol2inv(lm.i$qr$qr[p1, p1, drop=FALSE])
	    lhat[i,] <- t(x[i,]) %*% inv.Z %*% t(x) %*% diag(w.i)
        } else {
	    flag <- 1
	}
    }
    if (flag == 0) {
    	v1 <- sum(diag(lhat))
    	B1 <- t(diag(n)-lhat)%*%(diag(n)-lhat)
    	rss <- c(t(y)%*%B1%*%y)
    	sigma2.b <- rss / n
# NOTE 2* and sqrt() inserted for legibility
    	score <- 2*n*log(sqrt(sigma2.b)) + n*log(2*pi) + 
	    (n * (n + v1) / (n - 2 - v1))
    } else {
	score <- as.numeric(NA)
    }
    options(show.error.messages = TRUE)
    if (verbose) cat("Bandwidth:", q, "AIC:", score, "\n")
    score
}

gwr.cv.adapt.f <- function(q, y, x, coords, gweight, verbose=TRUE, longlat=FALSE)
{
    n <- NROW(x)
#    m <- NCOL(x)
    cv <- real(n)
    bw <- gw.adapt(dp=coords, fp=coords, quant=q, longlat=longlat)
    options(show.error.messages = FALSE)
    for (i in 1:n) {
        xx <- x[i, ]
	dxs <- spDistsN1(coords, coords[i,], longlat=longlat)
	if (!is.finite(dxs[i])) dxs[i] <- 0
	w.i <- gweight(dxs^2, bw[i])
        w.i[i] <- 0
	if (any(w.i < 0 | is.na(w.i)))
       		stop(paste("Invalid weights for i:", i))
        lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
        if(!inherits(lm.i, "try-error")) {
            b <- coefficients(lm.i)
            cv[i] <- y[i] - (t(b) %*% xx)
        }
    }
    score <- sqrt(sum(t(cv) %*% cv)/n)
    options(show.error.messages = TRUE)
    if (verbose) cat("Adaptive q:", q, "CV score:", score, "\n")
    score
}

