# Copyright 2001-2009 Roger Bivand and Danlin Yu
# 

gwr <- function(formula, data = list(), coords, bandwidth, 
	gweight=gwr.Gauss, adapt=NULL, hatmatrix=FALSE, fit.points, 
	longlat=FALSE, se.fit=FALSE, weights, cl=NULL, predictions=FALSE,
        fittedGWRobject=NULL, se.fit.CCT=TRUE) {
	this.call <- match.call()
	p4s <- as.character(NA)
	Polys <- NULL
	if (is(data, "SpatialPolygonsDataFrame")) 
		Polys <- as(data, "SpatialPolygons")
	if (is(data, "Spatial")) {
		if (!missing(coords))
		    warning("data is Spatial* object, ignoring coords argument")
		coords <- coordinates(data)
		p4s <- proj4string(data)
		data <- as(data, "data.frame")
	}

	if (missing(coords))
		stop("Observation coordinates have to be given")
	if (is.null(colnames(coords))) 
		colnames(coords) <- c("coord.x", "coord.y")

    	mf <- match.call(expand.dots = FALSE)
    	m <- match(c("formula", "data", "weights"), names(mf), 0)
    	mf <- mf[c(1, m)]
    	mf$drop.unused.levels <- TRUE
    	mf[[1]] <- as.name("model.frame")
    	mf <- eval(mf, parent.frame())
    	mt <- attr(mf, "terms")
	dp.n <- length(model.extract(mf, "response"))

    	weights <- as.vector(model.extract(mf, "weights"))
# set up default weights
    	if (!is.null(weights) && !is.numeric(weights)) 
        	stop("'weights' must be a numeric vector")
    	if (is.null(weights)) weights <- rep(as.numeric(1), dp.n)
    	if (any(is.na(weights))) stop("NAs in weights")
    	if (any(weights < 0)) stop("negative weights")
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)

	lm <- lm.wfit(x, y, w=weights)
	lm$x <- x
	lm$y <- y
	if (hatmatrix) se.fit <- TRUE
	if (hatmatrix) predictions <- TRUE

	if (missing(fit.points)) {
		fp.given <- FALSE
		fittedGWRobject <- NULL
		predictions <- TRUE
		fit.points <- coords
		colnames(fit.points) <- colnames(coords)
                predx <- x
	} else fp.given <- TRUE
	griddedObj <- FALSE
	if (is(fit.points, "Spatial")) {
		if (predictions) {
                    t1 <- try(slot(fit.points, "data"), silent=TRUE)
		    if (class(t1) == "try-error") 
			stop("No data slot in fit.points")
                    predx <- try(model.matrix(delete.response(mt), fit.points))
                    if (class(predx) == "try-error") 
			stop("missing RHS variable in fit.points")
                    if (ncol(predx) != ncol(x))
			stop("new data matrix columns mismatch")
                }
		Polys <- NULL
		if (is(fit.points, "SpatialPolygonsDataFrame")) {
			Polys <- Polygons(fit.points)
			fit.points <- coordinates(fit.points)
		} else {
			griddedObj <- gridded(fit.points)
			fit.points <- coordinates(fit.points)
		}
	} else {
        	if (predictions && fp.given)
		    stop("predictions not available for matrix fit points")
	}

	n <- NROW(fit.points)
	rownames(fit.points) <- NULL
	if (is.null(colnames(fit.points))) colnames(fit.points) <- c("x", "y")
        if (predictions) {
             if (nrow(predx) != nrow(fit.points))
		stop("new data matrix rows mismatch")
            fit.points <- cbind(fit.points, predx)
        }
#	if (is.null(fit.points)) fit.points <- coords
	m <- NCOL(x)
	if (NROW(x) != NROW(coords))
		stop("Input data and coordinates have different dimensions")
	if (missing(bandwidth) && is.null(adapt))
	    stop("Bandwidth must be given for non-adaptive weights")
	if (missing(bandwidth)) bandwidth <- NULL
	lhat <- NA
        yhat <- NULL
	if (!is.null(fittedGWRobject)) {
            yhat <- fittedGWRobject$SDF$pred
        }


	GWR_args <- list(fp.given=fp.given, hatmatrix=hatmatrix, 
	    longlat=longlat, bandwidth=bandwidth, adapt=adapt, se.fit=se.fit,
	    predictions=predictions, se.fit.CCT=se.fit.CCT)

	if (!is.null(cl) && length(cl) > 1 && fp.given && !hatmatrix) {
	    if (length(grep("cluster", class(cl))) > 0 && 
		.Platform$OS.type == "unix" && exists("splitIndices")) {
#		&& require(snow)) 
		l_fp <- lapply(splitIndices(nrow(fit.points), length(cl)), 
		    function(i) fit.points[i,])
		clusterEvalQ(cl, library(spgwr))

		clusterExport_l <- function(cl, list) {
                    gets <- function(n, v) {
                        assign(n, v, env = .GlobalEnv)
                        NULL
                    }
                    for (name in list) {
                        clusterCall(cl, gets, name, get(name))
                    }
		}

		clusterExport_l(cl, list("GWR_args", "coords", "gweight", "y",
		    "x", "weights", "yhat"))

		res <- parLapply(cl, l_fp, function(fp) .GWR_int(fit.points=fp,
		    coords=coords, gweight=gweight, y=y, x=x,
		    weights=weights, yhat=yhat, GWR_args=GWR_args))
		clusterEvalQ(cl, rm(list=c("GWR_args", "coords", "gweight", "y",
		    "x", "weights", "yhat")))
		df <- as.data.frame(do.call("rbind", 
		    lapply(res, function(x) x$df)))
		bw <- do.call("c", lapply(res, function(x) x$bw))
	        results <- NULL
	    } else stop("snow package not loaded")
	} else { # cl

	    df <- .GWR_int(fit.points=fit.points, coords=coords, 
		gweight=gweight, y=y, x=x, weights=weights, yhat=yhat,
		GWR_args=GWR_args)
	    if (!fp.given && hatmatrix) lhat <- df$lhat
	    bw <- df$bw
#	    df <- as.data.frame(df$df)
	    results <- NULL

	} # cl
	if (!fp.given && hatmatrix) {

	#This section calculates the effective degree of freedoms, edf;
	#the normalized residual sum of square, sigma2; the model
	#residual of squares, rss; and various version of AICs
	
	#As long as the hat matrix is obtained, many statistics can be 
	#calculated, such as the effective degree of freedom, AIC, etc.
	#Now calculate the effective degree of freedom of the residual
	#Reference: GWR book, page 55
	#obtain v1 and v2:
	
		v1 <- sum(diag(lhat))
		B2 <- t(lhat)%*%lhat
		v2 <- sum(diag(B2))
	
	#effective d.f. is n - 2*v1 + v2
	
		edf <- dp.n - 2*v1 + v2
	
	#Follow Leung et al. EPA 2000 page 15, the estimate of sigma square
	#can be obtained through rss and delta1 (which is actually edf)
	#Calculate rss:
	
		B1 <- t(diag(dp.n)-lhat)%*%(diag(dp.n)-lhat)
		rss <- c(t(y)%*%B1%*%y)
		delta1 <- sum (diag (B1))
		sigma2 <- rss/delta1 #line 77
	
	#Now the problem is, there are several version of AIC's calculation
	#formula, the GWR book's (page 61,96), Brunsdon's handouts, 
	#and the one from Hurvich, Simonoff and Tsai (1998, page 276)
	#I will implement all of them, and called them
	#AICb, from the book, AICh, from Brunsdon, and AICc, from Hurvich et al.
	
	#AICb <- n*log(sigma2) + n*log(2*3.14) + (n * (n + v1) / (n - 2 - v1))
	#AICh <- n*log(sigma2) + ((n + v1) / (n + 2 - v1))
	
	#To calculate AICc, there are several interal parameters
	#delta1 has already been calculated, detailed formula see
	#Hurvich et al. 1998, page 275, 276
	#B1 (above) and B2, delta2, nu1, nu2:
	

		odelta2 <- sum(diag(B1)^2)
# Patrick Zimmerman 090617
		delta2 <- sum(diag(B1 %*% B1))
		nu1 <- sum(diag(B2))
	#nu2 <- sum(diag(B2^)2)
	
	#AICc is from the formula in Hurvich et al. 1998 page 276
	#AICc1:
	
	#AICc <- n*log(sigma2) + n * ((delta1/delta2)*(n + nu1))/((delta1^2/delta2)-2)
	
	#One thing that I did not notice is that the sigma2 here I used is not
	#the same sigma2 used in the GWR book (detailed reference in page 96).
	#The sigma2 I used is from Leung et al (2000, p 15), calculated in line 77
	#The sigma2 in the GWR book is a maximum likelihood estimate
	#It should be: sigma2 <- rss/n instead of sigma2 <- rss/delta1
	#For this reason, a corrected AICb.b, AICh.b, AICc.b are therefore provided
	#followed by creating the sigma sqare used in the book, termed here sigma2.b
	#All the above unncessary calculation is then commented out.
	
		sigma2.b <- rss / dp.n
		AICb.b <- 2*dp.n*log(sqrt(sigma2.b)) + dp.n*log(2*pi) + 
			(dp.n * ((n + v1) / (dp.n - 2 - v1)))
# NOTE 2* and sqrt() inserted for legibility
		AICh.b <- 2*dp.n*log(sqrt(sigma2.b)) + dp.n*log(2*pi) + dp.n + v1
# added omitted n*log(2*pi) term in AICc.b
# bug resolved by Christian Salas 090418
		AICc.b <- 2*dp.n*log(sqrt(sigma2.b)) + dp.n*log(2*pi) + dp.n * 
			((delta1/delta2)*(dp.n + nu1))/((delta1^2/delta2)-2)
		results <- list(v1=v1, v2=v2, delta1=delta1, delta2=delta2, 
			sigma2=sigma2, sigma2.b=sigma2.b, AICb=AICb.b, 
			AICh=AICh.b, AICc=AICc.b, edf=edf, rss=rss, nu1=nu1,
                        odelta2=odelta2, n=dp.n)
	}
#	df <- data.frame(sum.w=sum.w, gwr.b, gwr.R2, gwr.se, gwr.e)
	if (!fp.given && is.null(fittedGWRobject)) {
	    localR2 <- numeric(n)	    
	    if (is.null(adapt)) {
		    bw <- bandwidth
		    bandwidthR2 <- rep(bandwidth, n)
	    } else {
		bandwidthR2 <- gw.adapt(dp=coords, fp=fit.points[,1:2], 
		    quant=adapt, longlat=longlat)
		bw <- bandwidthR2
	    }
	    if (any(bandwidth < 0)) stop("Invalid bandwidth")
	    for (i in 1:n) {
		dxs <- spDistsN1(coords, fit.points[i,1:2], 
		    longlat=GWR_args$longlat)
		if (any(!is.finite(dxs)))
			dxs[which(!is.finite(dxs))] <- 0
#		if (!is.finite(dxs[i])) dxs[i] <- 0
		w.i <- gweight(dxs^2, bandwidthR2[i])
		w.i <- w.i * weights
		if (any(w.i < 0 | is.na(w.i)))
        		stop(paste("Invalid weights for i:", i))
                RSS <- sum(w.i * (y - df$df[,"pred"])^2)
		yss <- sum(w.i * (y - weighted.mean(y, w.i))^2)
                localR2[i] <- 1 - (RSS/yss)
            }
            df$df <- cbind(df$df, localR2)
        }
        if (se.fit) {
            EDFS <- NULL
            normSigmaS <- NULL
            EDF <- NULL
            normSigma <- NULL
	    if (fp.given && !is.null(fittedGWRobject)) {
                if (fittedGWRobject$hatmatrix) {
		    EDF <- fittedGWRobject$results$edf
                    normSigma <- sqrt(fittedGWRobject$results$rss/EDF)
		    EDFS <- fittedGWRobject$results$n - 
                        fittedGWRobject$results$v1
                    normSigmaS <- sqrt(fittedGWRobject$results$rss/EDFS)
	        } 
	    }
	    if (!fp.given && hatmatrix) {
                EDFS <- results$n - results$v1
                normSigmaS <- sqrt(results$rss/EDFS)
                EDF <- results$edf
                normSigma <- sqrt(results$rss/EDF)
            }
            ses <- grep("_se", colnames(df$df))
            senms <- colnames(df$df)[ses]
            betase <- df$df[, ses]
            df$df[, ses] <- NA
            if (predictions) {
                pred.se <- df$df[, "pred.se"]
		df$df[, "pred.se"] <- NA
            }
            if (!is.null(EDF)) {
                betaseEDF <- normSigma * sqrt(betase)
                colnames(betaseEDF) <- paste(senms, "EDF", sep="_")
                df$df[, ses] <- normSigmaS * sqrt(betase)
                df$df <- cbind(df$df, betaseEDF)
                if (predictions) {
                    pred.se_EDF <- normSigma * sqrt(pred.se)
                    df$df[, "pred.se"] <- normSigmaS * sqrt(pred.se)
                    df$df <- cbind(df$df, pred.se_EDF)
                }
            } else {
                warning("standard errors set to NA, normalised RSS not available")
            }
        }

	df <- as.data.frame(df$df)
	if (predictions) fit.points <- fit.points[,1:2]
        row.names(fit.points) <- row.names(df)
	SDF <- SpatialPointsDataFrame(coords=fit.points, 
		data=df, proj4string=CRS(p4s))

	if (griddedObj) {
		gridded(SDF) <- TRUE
	} else {
		if (!is.null(Polys)) {
			df <- data.frame(SDF@data)
			rownames(df) <- sapply(slot(Polys, "polygons"),
                            function(i) slot(i, "ID"))
			SDF <- SpatialPolygonsDataFrame(Sr=Polys, data=df)
		}
	}
	z <- list(SDF=SDF, lhat=lhat, lm=lm, results=results, 
		bandwidth=bw, adapt=adapt, hatmatrix=hatmatrix, 
		gweight=deparse(substitute(gweight)), this.call=this.call)
	class(z) <- "gwr"
	invisible(z)
}


print.gwr <- function(x, ...) {
	if(class(x) != "gwr") stop("not a gwr object")
	cat("Call:\n")
	print(x$this.call)
	cat("Kernel function:", x$gweight, "\n")
	n <- nrow(coordinates(x$SDF))
	if (is.null(x$adapt)) cat("Fixed bandwidth:", x$bandwidth, "\n")
	else cat("Adaptive quantile: ", x$adapt, " (about ", 
		floor(x$adapt*n), " of ", n, ")\n", sep="")
	m <- length(x$lm$coefficients)
	cat("Summary of GWR coefficient estimates:\n")
	CM <- t(apply(as(x$SDF, "data.frame")[,(1+(1:m)), drop=FALSE], 2, summary))[,c(1:3,5,6)]
	if (is.null(dim(CM))) CM <- t(as.matrix(CM))
	CM <- cbind(CM, coefficients(x$lm))
	colnames(CM) <- c(colnames(CM)[1:5], "Global")
	printCoefmat(CM)
	if (x$hatmatrix) {
		cat("Number of data points:", n, "\n")
		cat("Effective number of parameters (residual: 2traceS - traceS'S):", 2*x$results$v1 -
		    x$results$v2, "\n")
		cat("Effective degrees of freedom (residual: 2traceS - traceS'S):", x$results$edf, "\n")
		cat("Sigma (residual: 2traceS - traceS'S):",
                    sqrt(x$results$rss/x$results$edf), "\n")
                cat("Effective number of parameters (model: traceS):",
                    x$results$v1, "\n")
		cat("Effective degrees of freedom (model: traceS):",
                    (x$results$n - x$results$v1), "\n")
		cat("Sigma (model: traceS):",
                    sqrt(x$results$rss/(x$results$n - x$results$v1)), "\n")
		cat("Sigma (ML):", sqrt(x$results$sigma2.b), "\n")
		cat("AICc (GWR p. 61, eq 2.33; p. 96, eq. 4.21):",
                    x$results$AICb, "\n")
		cat("AIC (GWR p. 96, eq. 4.22):", x$results$AICh, "\n")
		cat("Residual sum of squares:", x$results$rss, "\n")
	}
	invisible(x)
}

.GWR_int <- function(fit.points, coords, gweight, y, x, weights, yhat, 
	GWR_args) {
	    if (GWR_args$predictions) {
                predx <- fit.points[, -c(1,2)]
                fit.points <- fit.points[, c(1,2)]
            }
            
	    n <- nrow(fit.points)
 	    m <- NCOL(x)
            x1 <- matrix(1, nrow=nrow(x), ncol=1)
	    sum.w <- numeric(n)
            betas <- matrix(nrow=n, ncol=m)
	    colnames(betas) <- colnames(x)
            if(!GWR_args$fp.given) {
		gwr.e <- numeric(n)
	    } else {
		gwr.e <- NULL
	    }
	    if (GWR_args$se.fit) {
                betase <- matrix(nrow=n, ncol=m)
		colnames(betase) <- paste(colnames(x), "se", sep="_")
            } else {
                betase <- NULL
            }
            if (GWR_args$predictions) {
                pred <- numeric(n)
                if (GWR_args$se.fit) {
                    pred.se <- numeric(n)
                } else {
                    pred.se <- NULL
                }
            } else {
                pred <- NULL
		pred.se <- NULL
            }
	    if (is.null(yhat)) {
                localR2 <- NULL
            } else {
                localR2 <- numeric(n)
            }

	    if (!GWR_args$fp.given && GWR_args$hatmatrix) 
	        lhat <- matrix(nrow=n, ncol=n)
	    if (is.null(GWR_args$adapt)) {
		    bw <- GWR_args$bandwidth
		    bandwidth <- rep(GWR_args$bandwidth, n)
	    } else {
		bandwidth <- gw.adapt(dp=coords, fp=fit.points, 
		    quant=GWR_args$adapt, longlat=GWR_args$longlat)
		bw <- bandwidth
	    }
	    if (any(bandwidth < 0)) stop("Invalid bandwidth")
	    for (i in 1:n) {
		dxs <- spDistsN1(coords, fit.points[i,], 
		    longlat=GWR_args$longlat)
		if (any(!is.finite(dxs)))
			dxs[which(!is.finite(dxs))] <- 0
#		if (!is.finite(dxs[i])) dxs[i] <- 0
		w.i <- gweight(dxs^2, bandwidth[i])
		w.i <- w.i * weights
		if (any(w.i < 0 | is.na(w.i)))
        		stop(paste("Invalid weights for i:", i))
		lm.i <- lm.wfit(x, y, w.i)
		sum.w[i] <- sum(w.i)
		betas[i,] <- coefficients(lm.i)
		ei <- residuals(lm.i)
# prediction fitted values at fit point
                if (GWR_args$predictions) {
                    pred[i] <- sum(predx[i,] * betas[i,])
                }
# use of diag(w.i) dropped to avoid forming n by n matrix
# bug report: Ilka Afonso Reis, July 2005

# local R-squared as explained by Tomoki Nakaya, 090624 (in GWR4)
# differs from local weighted regression R-squared

		if (!is.null(yhat)) {
		    RSS <- sum(w.i * (y - yhat)^2)
                    yss <- sum(w.i * (y - weighted.mean(y, w.i))^2)
		    localR2[i] <- 1 - (RSS/yss)
		}

	        if (GWR_args$se.fit) {
		    p <- lm.i$rank
		    p1 <- 1:p
		    inv.Z <- chol2inv(lm.i$qr$qr[p1, p1, drop=FALSE])
# p. 55 CC definition 
                    if (GWR_args$se.fit.CCT) {
                        C <- inv.Z %*% t(x) %*% diag(w.i)
                        CC <- C %*% t(C)
# only return coefficient covariance matrix diagonal raw values
# for post-processing
                        betase[i,] <- diag(CC)
                    } else {
                        betase[i,] <- diag(inv.Z)
                    }
# prediction "standard errors"
# only return raw values for post-processing
                    if (GWR_args$predictions) {
                        if (GWR_args$se.fit.CCT) {
                            pred.se[i] <- t(predx[i,]) %*% CC %*% predx[i,]
                        } else {
                            pred.se[i] <- t(predx[i,]) %*% inv.Z %*% predx[i,]
                        }
		    }
		}
# assigning residual bug Torleif Markussen Lunde 090529
		if (!GWR_args$fp.given) gwr.e[i] <- ei[i]

		if (!GWR_args$fp.given && GWR_args$hatmatrix) 
			lhat[i,] <- t(x[i,]) %*% inv.Z %*% t(x) %*% diag(w.i)
	    }
	    df <- cbind(sum.w, betas, betase, gwr.e, pred, pred.se, localR2)
	    if (!GWR_args$fp.given && GWR_args$hatmatrix) 
		return(list(df=df, lhat=lhat, bw=bw))
	    else return(list(df=df, bw=bw))
} # GWR_int


