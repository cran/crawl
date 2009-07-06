\name{crwSamplePar}

\alias{crwSamplePar}

% -------------------------------------------------------------------------

\title{
   Create a weighted importance sample for posterior predictive track simulation.
}


\description{

  The \code{crwSamplePar} function uses a fitted model object from
  \code{crwMLE} and a set of prediction times to construct a list from which
  \code{\link{crwPostIS}} will draw a sample from either the posterior
  distribution of the state vectors conditional on fitted parameters or a full
  posterior draw from an importance sample of the parameters.

}


\usage{
crwSamplePar(object.sim, size=1000, df=Inf, scale=1)

}


\arguments{
  \item{object.sim}{A simulation object from \code{\link{crwSimulator}}.}

  
  \item{size}{Size of the parameter importance sample}
  
  \item{df}{Degrees of freedom for the t approximation to the parameter posterior}
  
  \item{scale}{Scale multiplier for the covariance matrix of the t approximation}

}


\value{

  List with the following elements:

  \item{x}{Longitude coordinate with NA at prediction times}

  \item{y}{Similar to above for latitude}

  \item{locType}{Indicates prediction types with a "p" or observation times with
                 an "o"}
  \item{P1.y}{Initial state covariance for latitude}
  
  \item{P1.x}{Initial state covariance for longitude}
  
  \item{a1.y}{Initial latitude state}
  
  \item{a1.x}{Initial longitude state}
  
  \item{n.errX}{number of longitude error model parameters}
  
  \item{n.errY}{number of latitude error model parameters}
  
  \item{delta}{vector of time differences}
  
  \item{driftMod}{Logical. indicates random drift model}
  
  \item{stopMod}{Logical. Indicated stop model fitted}
  
  \item{stop.mf}{stop model design matrix}
  
  \item{err.mfX}{Longitude error model design matrix}
  
  \item{err.mfY}{Latitude error model design matrix}
  
  \item{mov.mf}{Movement model design matrix}
  
  \item{fixPar}{Fixed values for parameters in model fitting}
  
  \item{Cmat}{Covaraince matrix for parameter sampling distribution}
  
  \item{Lmat}{Cholesky decomposition of Cmat}
  
  \item{par}{fitted parameter values}
  
  \item{N}{Total number of locations}
  
  \item{loglik}{log likelihood of the fitted model}
  
  \item{Time}{vector of observation times}
  
  \item{coord}{names of coordinate vectors in original data}
  
  \item{Time.name}{Name of the observation times vector in the original data}
  
  \item{thetaSampList}{A list containing a data frame of parameter vectors and their associated
                       probabilities for a resample}
  

}


\details{

  The crwSimulator function produces a list and preprocesses the necessary
  components for repeated track simulation from a fitted CTCRW model from
  \code{\link{crwMLE}}. If df=Inf (default) then a multivariate normal distribution
  is used to approximate the parameter posterior.  

}


\seealso{ See \code{\link{northernFurSeal}} for example. }

% -------------------------------------------------------------------------

\author{Devin S. Johnson}
