\name{crwPredict}

\alias{crwPredict}

% -------------------------------------------------------------------------

\title{ Predict animal locations and velocities using a fitted CTCRW
  model }


\description{

  The \code{crwPredict} function uses a fitted model object from
  \code{crwMLE} to predict animal locations (with estimated uncertainty)
  at times in the original data set and supplimented by times in
  \code{predTime}. If \code{speedEst} is set to \code{TRUE}, then animal
  log-speed is also estimated.

}


\usage{
crwPredict(object.crwFit, predTime=NULL, speedEst=FALSE, flat=FALSE)

}


\arguments{
  \item{object.crwFit}{A model object from \code{\link{crwMLE}}.}

  \item{predTime}{vector of additional prediction times.}

  \item{speedEst}{logical. Estimate animal speed or not.}

  \item{flat}{logical. Should the result be returned as a flat
    data.frame.}

}


\value{

  List with the following elements:

  \item{originalData}{A data.frame with is \code{data} merged with
    \code{predTime}.}

  \item{alpha.hat.y}{A data.frame with predicted state values for each
    time. First column in latitude location (\code{mu.y}), second in
    velocity (\code{nu.y} or \code{theta.y} for drift models), and third
    is drift velocity (\code{gamma.y} if specified).}

  \item{alpha.hat.x}{longitude state predictions.}

  \item{Var.hat.y}{array where \code{Var.hat.y[,,i]} is the prediction
    covariance matrix for \code{alpha.hat.y[,i]}.}

  \item{Var.hat.x}{array or covariance matrices for \code{alpha.hat.x}.}

  \item{speed}{(If \code{speedEst=TRUE}) Gives log speed estimates for
  each time and standard errors based on delta method. If coordinates
  are polar, units are meters/unit \code{Time}, else, units are those
  specified by the coordinates.}

If \code{flat} is set to \code{TRUE} then a data set is returned with
the columns of the original data plus the state estimates, standard
errors (se) and the speed estimates.

}


\details{

  The requirements for \code{data} are the same as those for fitting the
  model in \code{\link{crwMLE}}.

}


\seealso{ See \code{\link{northernFurSeal}} for example. }

% -------------------------------------------------------------------------

\author{Devin S. Johnson}
