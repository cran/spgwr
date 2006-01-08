require(sp)
.First.lib <- function(lib, pkg) {
	library.dynam("spgwr", pkg, lib)
}
