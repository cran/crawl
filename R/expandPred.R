"expandPred" <- function(origTime, x, predTime)
{
    newx <- merge(data.frame(Time=origTime, x), data.frame(Time=predTime),
                  by=c('Time'), all=TRUE)
    out <- apply(as.matrix(newx[, -1]), 2, function(vec, Time) {
        approx(Time, vec, Time, 'const')$y
    }, Time=newx[, 1])
    return(out[!duplicated(newx$Time), ])
}
