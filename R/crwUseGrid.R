# TODO: Add comment
# 
# Author: johnsond@afsc.noaa.gov
###############################################################################

crwUseGrid <- function(object, grid, rm.zeros=FALSE){
	if(!inherits(object, "crwPredict")) stop("The argument 'object' must be of class crwPredict!")
	if(!inherits(grid,"GridTopology")) stop("Argument 'grid' must be a 'GridTopology' object!")
	firstCell <- as.vector(grid@cellcentre.offset)
	eps <- as.vector(grid@cellsize)
	dims <- as.vector(grid@cells.dim)
	x.center <- seq(from=firstCell[1], by=eps[1], length=dims[1])
	y.center <- seq(from=firstCell[2], by=eps[2], length=dims[2])
	xy.centers <- list(x=x.center, y=y.center)	
	x <- sort(xy.centers$x)
	y <- sort(xy.centers$y)
	x.pts <- object$mu.x[object$locType=="p"]
	y.pts <- object$mu.y[object$locType=="p"]
	eps.x <- min(diff(x))
	eps.y <- min(diff(y))
	x.cut <- c(x[1]-eps.x/2, x+eps.x/2)
	y.cut <- c(y[1]-eps.y/2, y+eps.y/2)
	mat <- matrix(0,length(y),length(x))
	x.int <- findInterval(x.pts, x.cut)
	y.int <- findInterval(y.pts, y.cut)
	for(i in 1:length(x.pts)){
		mat[length(y)-y.int[i]+1, x.int[i]] <- mat[length(y)-y.int[i]+1, x.int[i]] + 1
	}
	out <- SpatialGridDataFrame(grid, data.frame(use = as.vector(t(mat))))
	if(rm.zeros) out$use <- ifelse(out$use==0, NA, out$use)
	return(out)
}



