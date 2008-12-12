"crwPredict" <- function(object.crwFit, data, predTime=NULL,
                         speedEst=FALSE, flat=FALSE)
{
    ## Model definition/parameters ##
    driftMod <- object.crwFit$random.drift
    stopMod <- !is.null(object.crwFit$stop.model)
    mov.mf <- object.crwFit$mov.mf
    stop.mf <- object.crwFit$stop.mf
    err.mfX <- object.crwFit$err.mfX
    err.mfY <- object.crwFit$err.mfY
    par <- object.crwFit$par
    n.errX <- object.crwFit$n.errX
    n.errY <- object.crwFit$n.errY
    n.mov <- object.crwFit$n.mov
    tn <- object.crwFit$Time.name

    ## Data setup ##
    if (!is.null(predTime)) {
        origTime <- data[,tn]
        if (is.null(data$locType)) data$locType <- "o"
        predData <- data.frame(predTime, "p")
        names(predData) <- c(tn, "locType")
        data <- merge(data, predData,
                      by=c(tn, "locType"), all=TRUE)
        dups <- duplicated(data[,tn])
        data <- data[!dups, ]
        mov.mf <- as.matrix(expandPred(origTime, mov.mf, predTime))
        if (stopMod) stop.mf <- as.matrix(expandPred(origTime, stop.mf, predTime))
        if (!is.null(err.mfX)) err.mfX <- as.matrix(expandPred(origTime, err.mfX, predTime))
        if (!is.null(err.mfY)) err.mfY <- as.matrix(expandPred(origTime, err.mfY, predTime))
    }
    if (object.crwFit$polar.coord) {
        lonAdjVals <- cos(round(approx(data[, tn], data[, object.crwFit$coord[2]],
                                       data[, tn])$y, 0) * pi / 180)
    } else lonAdjVals <- rep(1, nrow(data))
    delta <- c(diff(data[, tn]), 1)
    a1.x <- object.crwFit$initial.state$a1.x
    P1.x <- object.crwFit$initial.state$P1.x
    a1.y <- object.crwFit$initial.state$a1.y
    P1.y <- object.crwFit$initial.state$P1.y
    y <- data[, object.crwFit$coord[2]]
    x <- data[, object.crwFit$coord[1]]
    loctype <- ifelse(is.na(x) | is.na(y), 1, 0)
    y <- ifelse(loctype == 1, 9999, y)
    x <- ifelse(loctype == 1, 9999, x)

    if (!is.null(err.mfX)) {
        theta.errX <- par[1:n.errX]
        tau2x <- exp(2 * err.mfX %*% theta.errX)
    } else tau2x <- 0.0
    if (!is.null(err.mfY)) {
        theta.errY <- par[(n.errX + 1):(n.errX + n.errY)]
        tau2y <- exp(2 * err.mfY %*% theta.errY)
    } else tau2y <- tau2x
    theta.mov <- par[(n.errX + n.errY + 1):(n.errX + n.errY + 2 * n.mov)]
    sig2 <- exp(2 * (mov.mf %*% theta.mov[1:n.mov]))
    b <- exp(2 * (mov.mf %*% theta.mov[(n.mov + 1):(2 * n.mov)]))
    stay <- rep(0, nrow(data))
    if (stopMod) {
        stop.mf <- object.crwFit$stop.mf
        theta.stop <- par[(n.errX + n.errY + 2 * n.mov + 1)]
        b[stop.mf != 1] <- b[stop.mf != 1] / ((1 - stop.mf[stop.mf != 1]) ^
                               exp(theta.stop))
        stay[stop.mf == 1] <- 1
    }
    if (driftMod) {
        theta.drift <- par[(n.errX + n.errY + 2 * n.mov + 1):
                           (n.errX + n.errY + 2 * n.mov + 2)]
        b.drift <- b / exp(theta.drift[2])
        sig2.drift <- exp(2 * theta.drift[1])
        call.predict <- "crwDrift_predict"
    } else {
        b.drift <- sig2.drift <- 0.0
        call.predict <- "crw_predict"
    }
    N <- nrow(data)
    out <- .Fortran(call.predict,
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
                    lonadj=as.double(lonAdjVals),
                    N=as.integer(N),
                    lly=as.double(0),
                    llx=as.double(0),
                    predy=as.double(matrix(0, N, 2 + driftMod)),
                    predx=as.double(matrix(0, N, 2 + driftMod)),
                    vary=as.double(array(0, c(2 + driftMod, 2 + driftMod, N))),
                    varx=as.double(array(0, c(2 + driftMod, 2 + driftMod, N))),
                    package="crawl")
    predy <- data.frame(matrix(out$predy, N, 2 + driftMod))
    if (driftMod) {
        names(predy) <- c("mu.y", "theta.y", "gamma.y")
    } else names(predy) <- c("mu.y", "nu.y")
    predx <- data.frame(matrix(out$predx, N, 2 + driftMod))
    if (driftMod) {
        names(predx) <- c("mu.x", "theta.x", "gamma.x")
    } else names(predx) <- c("mu.x", "nu.x")
    vary <- array(out$vary, c(2 + driftMod, 2 + driftMod, N))
    varx <- array(out$varx, c(2 + driftMod, 2 + driftMod, N))
    if (speedEst) {
        log.speed <- logSpeed(predx, predy, varx, vary, object.crwFit$polar.coord)
    } else log.speed <- NULL
    out <- list(originalData=data, alpha.hat.y=predy, alpha.hat.x=predx,
                V.hat.y=vary, V.hat.x=varx, speed=log.speed, loglik=out$lly + out$llx)
    attr(out, "coord") <- c(x=object.crwFit$coord[1], y=object.crwFit$coord[2])
    attr(out, "random.drift") <- driftMod
    attr(out, "stop.model") <- object.crwFit$stopMod
    attr(out, "polar.coord") <- object.crwFit$polar.coord
    attr(out, "Time.name") <- tn
    if (flat) {
        out <- as.flat(out)
    } else {
        class(out) <- c("crwPredict", "list")
        attr(out, "flat") <- FALSE
    }
    return(out)
}
