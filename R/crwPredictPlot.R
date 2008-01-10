"crwPredictPlot" <- function(object, plotType="ll")
{
    y.c <- attr(object, "coord")['y']
    x.c <- attr(object, "coord")['x']
    if (!attr(object, "flat")) {
        mu.xUp <- object$alpha.hat.x[, 1] + 1.96 * sqrt(object$V.hat.x[1, 1, ])
        mu.xLo <- object$alpha.hat.x[, 1] - 1.96 * sqrt(object$V.hat.x[1, 1, ])
        mu.yUp <- object$alpha.hat.y[, 1] + 1.96 * sqrt(object$V.hat.y[1, 1, ])
        mu.yLo <- object$alpha.hat.y[, 1] - 1.96 * sqrt(object$V.hat.y[1, 1, ])
        xvals <- object$originalData[, x.c]
        yvals <- object$originalData[, y.c]
        mu.x <- object$alpha.hat.x[, 1]
        mu.y <- object$alpha.hat.y[, 1]
        Time <- object$originalData[, attr(object, "Time.name")]
    } else {
        mu.xUp <- object$mu.x + 1.96 * object$se.mu.x
        mu.xLo <- object$mu.x - 1.96 * object$se.mu.x
        mu.yUp <- object$mu.y + 1.96 * object$se.mu.y
        mu.yLo <- object$mu.y - 1.96 * object$se.mu.y
        xvals <- object[, x.c]
        yvals <- object[, y.c]
        mu.x <- object$mu.x
        mu.y <- object$mu.y
        Time <- object[, attr(object, "Time.name")]
    }

    mu.y.mx <- max(pmax(mu.yUp, yvals, na.rm=TRUE))
    mu.y.mn <- min(pmin(mu.yLo, yvals, na.rm=TRUE))
    mu.x.mx <- max(pmax(mu.xUp, xvals, na.rm=TRUE))
    mu.x.mn <- min(pmin(mu.xLo, xvals, na.rm=TRUE))
    y.ylims <- c(mu.y.mn, mu.y.mx)
    x.ylims <- c(mu.x.mn, mu.x.mx)

    switch(plotType,
           map = {
               plot(xvals, yvals, pch=16, col="blue", xlab=x.c, ylab=y.c,
                    xlim=x.ylims, ylim=y.ylims, cex=0.5)
               lines(mu.x, mu.y, col="red")},
           ll = {layout(matrix(1:2, ncol=1))
                 plot(Time, xvals, pch=16, col="blue", xlab="time", ylab=x.c,
                      ylim=x.ylims, cex=0.5)
                 lines(Time, mu.x, col="red")
                 lines(Time, mu.xUp, col="green", pch=16, cex=0.2)
                 lines(Time, mu.xLo, col="green", pch=16, cex=0.2)
                 plot(Time, yvals, pch=16, col="blue", xlab="time", ylab=y.c,
                      ylim=y.ylims, cex=0.5)
                 lines(Time, mu.y, col="red")
                 lines(Time, mu.yUp, col="green", pch=16, cex=0.2)
                 lines(Time,mu.yLo, col='green', pch=16, cex=0.2)})
}