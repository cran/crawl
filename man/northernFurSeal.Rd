\name{northernFurSeal}
\docType{data}

\alias{northernFurSeal}

% -------------------------------------------------------------------------

\title{Northern fur seal pup telemetry data set}


\description{

  Northern fur seal pup relocation data set used in Johnson et
  al. (2008)

}


\usage{

data(northernFurSeal)

}


\format{
  A data frame with 795 observations on the following 4 variables:

  \describe{
    \item{\code{Time}}{a numeric vector.}

    \item{\code{Argos_loc_class}}{a factor with levels \code{0} \code{1}
      \code{2} \code{3} \code{A}.}

    \item{\code{latitude}}{a numeric vector.}

    \item{\code{longitude}}{a numeric vector.}

  }

}


\source{
  Alska Ecosystems Program
  National Marine Mammal Laboratory
  Alaska Fisheries Science Center
  National Marine Fisheries Service, NOAA
  7600 Sand Point Way NE
  Seattle, WA 98115
}


\references{

  Johnson, D., J. London, M. -A. Lea, and J. Durban (2008)
  Continuous-time random walk model for animal telemetry data. Ecology
  89:1208-1215.

}

% -------------------------------------------------------------------------

\examples{

data(northernFurSeal)

argosClasses <- c("3", "2", "1", "0", "A", "B")
ArgosMultFactors <- data.frame(Argos_loc_class=argosClasses,
                               errX=log(c(1, 1.5, 4, 14, 5.21, 20.78)),
                               errY=log(c(1, 1.5, 4, 14, 11.08, 31.03)))
nfsNew <- merge(northernFurSeal, ArgosMultFactors,
                by=c("Argos_loc_class"), all.x=TRUE)
nfsNew <- nfsNew[order(nfsNew$Time), ]

# State starting values
initial.drift <- list(a1.x=c(189.686, 0, 0), a1.y=c(57.145, 0, 0),
                      P1.x=diag(c(0, 0.001, 0.001)),
                      P1.y=diag(c(0, 0.001, 0.001)))

##Fit random drift model
fit <- crwMLE(mov.model=~1, err.model=list(x=~errX, y=~errY), drift.model=TRUE,
              data=nfsNew, coord=c("longitude", "latitude"),
              Time.name="Time", initial.state=initial.drift,
              fixPar=c(NA, 1, NA, 1, NA, NA, NA,NA),
              control=list(maxit=2000,trace=1, REPORT=10))

##Make hourly location predictions
predTime <- seq(ceiling(min(nfsNew$Time)), floor(max(nfsNew$Time)), 1)
predObj <- crwPredict(object.crwFit=fit, predTime=predTime, speedEst=TRUE,
                      flat=FALSE)
str(predObj)
predObj <- as.flat(predObj)
head(predObj)
crwPredictPlot(predObj)

}


\keyword{datasets}
