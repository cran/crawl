\name{expandPred}

\alias{expandPred}

% -------------------------------------------------------------------------

\title{Expand a time indexed data set with additional prediction times}


\description{

  Expands a covariate data frame (or vector) that has a separate time
  index by inserting prediction times and duplicating the covariate
  values for all prediction time between subsequent data times.

}


\usage{

expandPred(origTime, x, predTime)

}


\arguments{

  \item{origTime}{Original time index for the data to be expanded.}

  \item{x}{Data to be expanded.}

  \item{predTime}{prediction times to expand data.}

}


\value{ data.frame expanded by \code{predTime}. }

% -------------------------------------------------------------------------

\examples{

origTime <- c(1:10)
x <- cbind(rnorm(10), c(21:30))
predTime <- seq(1,10, by=0.25)
expandPred(origTime, x, predTime)

}


\author{Devin S. Johnson}
