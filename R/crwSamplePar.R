crwSamplePar <- function(object.sim, size=1000, df=Inf, scale=1)
{
   if(!inherits(object.sim, 'crwSimulator'))
      stop("Argument needs to be of class 'crwSimulator'\nUse 'crwSimulator( )' to create")
   fixPar <- object.sim$fixPar
   Cmat <- object.sim$Cmat[is.na(fixPar),is.na(fixPar)]
   se <- sqrt(diag(Cmat))
   err.mfX <- object.sim$err.mfX
   err.mfY <- object.sim$err.mfY
   parMLE <- object.sim$par
   stopMod <- object.sim$stopMod
   driftMod <- object.sim$driftMod
   stop.mf <- object.sim$stop.mf
   mov.mf <- object.sim$mov.mf
   n.errX <- object.sim$n.errX
   n.errY <- object.sim$n.errY
   n.mov <- object.sim$n.mov
   N <- object.sim$N
   lower <- object.sim$lower
   upper <- object.sim$upper 
   delta <- object.sim$delta
   thetaMat <- matrix(NA, size, length(fixPar)+3)
   for(i in 1:(size-1)){
   	 par <- parMLE
     eInd <- is.na(fixPar)
   	 eps <- rmvtt(mu=rep(0,sum(eInd)), Sigma=scale*Cmat, df=df, lower-par[eInd], upper-par[eInd])
     par[eInd] <- parMLE[eInd] + eps
     if(df==Inf) dens <- dmvnorm(eps, sigma=scale*Cmat, log=TRUE) - dmvnorm(0.0*eps, sigma=scale*Cmat, log=TRUE)
     else dens <- dmvt(eps, sigma=scale*Cmat, df=df, log=TRUE) - dmvt(0.0*eps, sigma=scale*Cmat, df=df, log=TRUE)
   ###
   ### Process parameters for Fortran
   ###
   if (!is.null(err.mfX)) {
      theta.errX <- par[1:n.errX]
      tau2x <- exp(2 * err.mfX %*% theta.errX)
   } else tau2x <- rep(0.0, N)
   if (!is.null(err.mfY)) {
      theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
      tau2y <- exp(2 * err.mfY %*% theta.errY)
   } else tau2y <- tau2x
   theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
   sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
   b <- exp(mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)])
   stay <- rep(0, N)
   if (stopMod) {
      stop.mf <- stop.mf
      theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
      b <- b / ((1 - stop.mf) ^ exp(theta.stop))
      stay <- ifelse(b==Inf, 1, 0)
      b <- ifelse(b==Inf, 9999, b) 
   }
   if (driftMod) {
      theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                                    (n.errX + n.errY + 2 * n.mov + 2)]
      b.drift <- exp(log(b) - log(1+exp(theta.drift[2])))
      sig2.drift <- exp(log(sig2) + 2 * theta.drift[1])
      call.lik <- "crwdriftn2ll"
   } else {
      b.drift <- sig2.drift <- 0.0
      call.lik <- "crwn2ll"
   }
   movMats <- getQT(sig2, b, sig2.drift, b.drift, delta, driftMod)
     out <- .Fortran(call.lik,
                    tau2y=as.double(tau2y),
                    tau2x=as.double(tau2x),
                    Qmat=as.double(movMats$Qmat),
                    Tmat=as.double(movMats$Tmat),
                    x=as.double(object.sim$x),
                    y=as.double(object.sim$y),
                    loctype=as.integer(object.sim$loctype),
                    stay=as.integer(stay),
                    ay=as.double(object.sim$a1.y),
                    ax=as.double(object.sim$a1.x),
                    Py=as.double(object.sim$P1.y),
                    Px=as.double(object.sim$P1.x),
                    lonadj=as.double(object.sim$lonAdj),
                    N=as.integer(N),
                    lly=as.double(0),
                    llx=as.double(0),
                    package="crawl")
     thetaMat[i,] <- c(out$lly+out$llx - dens, out$lly+out$llx, dens, par)
  }
thetaMat[size,] <- c(object.sim$loglik, object.sim$loglik, 0, object.sim$par)
thetaMat[,1] <- exp(thetaMat[,1]-max(thetaMat[,1]))/sum(exp(thetaMat[,1]-max(thetaMat[,1])))
colnames(thetaMat) <- c("w", "lik", "prop.lik", object.sim$nms)
attr(thetaMat,"effSamp") <- size/(1+(sd(thetaMat[,"w"])/mean(thetaMat[,"w"]))^2) 
if(is.null(object.sim$thetaSampList)) object.sim$thetaSampList <- list(thetaMat)
else object.sim$thetaSampList <- append(object.sim$thetaSampList, list(thetaMat))
return(object.sim)	
}
     