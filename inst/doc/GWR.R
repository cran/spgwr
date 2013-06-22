### R code from vignette source 'GWR.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: GWR.Rnw:38-43
###################################################
owidth <- getOption("width")
options("width"=90)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 3: afigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")


###################################################
### code chunk number 4: bfigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")


###################################################
### code chunk number 5: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 6: zfig (eval = FALSE)
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 7: zfigl (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 8: GWR.Rnw:103-109
###################################################
library(maptools)
library(spdep)
owd <- getwd()
setwd(system.file("etc/shapes", package="spdep"))
NY8 <- readShapeSpatial("NY8_utm18")
setwd(owd)


###################################################
### code chunk number 9: GWR.Rnw:115-117
###################################################
load(system.file("backstore/nyGWR.RData", package="spgwr"))
#save(bwG, gwrG, gbwG, ggwrG, file="nyGWR.RData")


###################################################
### code chunk number 10: GWR.Rnw:120-121
###################################################
library(spgwr)


###################################################
### code chunk number 11: GWR.Rnw:123-125 (eval = FALSE)
###################################################
## bwG <- gwr.sel(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, gweight=gwr.Gauss, verbose=FALSE)
## gwrG <- gwr(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, bandwidth=bwG, gweight=gwr.Gauss, hatmatrix=TRUE)


###################################################
### code chunk number 12: GWR.Rnw:127-128
###################################################
gwrG


###################################################
### code chunk number 13: GWR.Rnw:164-166 (eval = FALSE)
###################################################
## gbwG <- ggwr.sel(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", gweight=gwr.Gauss, verbose=FALSE)
## ggwrG <- ggwr(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8, family="poisson", bandwidth=gbwG, gweight=gwr.Gauss)


###################################################
### code chunk number 14: GWR.Rnw:168-169 (eval = FALSE)
###################################################
## ggwrG


###################################################
### code chunk number 15: GWR.Rnw:171-175
###################################################
xx <- capture.output(print(ggwrG))
cat(xx[1], sep="\n")
cat(strwrap(paste(xx[2:3], collapse=" "), width=70, exdent=4), sep="\n")
cat(xx[4:11], sep="\n")


###################################################
### code chunk number 16: GWR.Rnw:190-193
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
spplot(ggwrG$SDF, "PEXPOSURE", col.regions=grey.colors(7, 0.95, 0.55, 2.2), cuts=6)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 17: GWR.Rnw:198-201
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
pairs(as(ggwrG$SDF, "data.frame")[,2:5])
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 18: GWR.Rnw:209-211
###################################################
options("width"=owidth)
options("warn"=ow)


