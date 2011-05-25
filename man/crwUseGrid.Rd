\name{crwUseGrid}
\alias{crwUseGrid}
\title{
Compute a spatial use grid from a crawl prediction
}
\description{
This function take a crwPredict object and a spatial GridTopology object from the 'sp' package 
and outputs the number of predicted locations in each grid cell
}
\usage{
crwUseGrid(object, grid, rm.zeros = FALSE)
}

\arguments{
  \item{object}{
A crwPredict object output from \code{\link{crwPredict}}. The prediction object MUST have been created 
using a set of prediction times (i.e., is.null(predTime)=FALSE)
}
  \item{grid}{
A \code{GridTopology} object from the 'sp' package
}
  \item{rm.zeros}{
Logical. If set to true, the zeros in the use grid are set to NA for better plots. Use grids 
can have an overwhelming number of zeros.
}
}

\value{
A \code{SpatialGridDataFrame} with data column 'use' which gives the count of predicted (see the 'sp' package)
locations within the grid cell
}

\author{
Devin S. Johnson <devin.johnson@noaa.gov>
}
