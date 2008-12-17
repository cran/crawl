"crwN2ll" <- function(theta, fixPar, y, x, loctype, delta, a1.y, a1.x,
                      P1.x, P1.y, lonAdj, mov.mf, err.mfX, err.mfY, stop.mf,
                      n.errX, n.errY, n.mov, stopMod, driftMod)
{
    N <- length(y)
    fixPar[is.na(fixPar)] <- theta
    if (!is.null(err.mfX)) {
        theta.errX <- fixPar[1:n.errX]
        tau2x <- exp(2 * err.mfX %*% theta.errX)
    } else tau2x <- 0.0
    if (!is.null(err.mfY)) {
        theta.errY <- fixPar[(n.errX + 1):(n.errX + n.errY)]
        tau2y <- exp(2 * err.mfY %*% theta.errY)
    } else tau2y <- tau2x
    theta.mov <- fixPar[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
    sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
    b <- exp(2 * (mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)]))
    stay <- rep(0, N)
    if (stopMod) {
        theta.stop <- fixPar[(n.errX + n.errY + 2 * n.mov + 1)]
        b[stop.mf != 1] <- b[stop.mf != 1] / ((1 - stop.mf[stop.mf != 1]) ^
                               exp(theta.stop))
        stay[stop.mf == 1] <- 1
    }
    if (driftMod) {
        theta.drift <- fixPar[(n.errX + n.errY + 2 * n.mov + 1):
                              (n.errX + n.errY + 2 * n.mov + 2)]
        b.drift <- b / exp(theta.drift[2])
        sig2.drift <- exp(2 * theta.drift[1])
        call.lik <- "crwdriftn2ll"
    } else {
        b.drift <- sig2.drift <- 0.0
        call.lik <- "crwn2ll"
    }
    out <- .Fortran(call.lik,
                    tau2y=as.double(tau2y),
                    tau2x=as.double(tau2x),
                    sig2=as.double(sig2),
                    b=as.double(b),
                    bd=as.double(b.drift),
                    sig2d=as.double(sig2.drift),
                    delta=as.double(delta),
                    x=as.double(x),
                    y=as.double(y),
                    loctype=as.integer(loctype),
                    stay=as.integer(stay),
                    ay=as.double(a1.y),
                    ax=as.double(a1.x),
                    Py=as.double(P1.y),
                    Px=as.double(P1.x),
                    lonadj=as.double(lonAdj),
                    N=as.integer(N),
                    lly=as.double(0),
                    llx=as.double(0),
                    package="crawl")
    return(-2 * (out$lly + out$llx))
}
