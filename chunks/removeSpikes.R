## remove spikes
## Criterion: deviation more than 5 times running-average standard deviation
removeSpikes <- function (v, sdLimit=4) {
  mean <- rollapply(v, width=99, FUN=mean, fill=NA)
  sd <- rollapply(v, width=99, FUN=sd, fill=NA)
  # plotWAC(data.frame(DataHR$Time, (DataHR$QCF-mean)/sd))
  ix <- which(abs(v - mean) / sd > sdLimit)
  v[ix] <- NA
  v <- SmoothInterp (v, .Length=0)
  return(v)
}
