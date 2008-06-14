require(sp)
.First.lib <- function(lib, pkg) {
	library.dynam("spgwr", pkg, lib)
	packageStartupMessage(paste(
	"NOTE: default kernel and CV criteria changed", 
	"see help pages for details", sep="\n"), appendLF = FALSE)
}
