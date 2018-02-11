
# merge Matt's file (name like ARISTO2017rf04_LAMS_TAS_Update.nc) with standard to get
# rf04_MH.nc

Project <- 'ARISTO2017'
source ('./getNetCDFMH.R')
transferAttributes <- function (d, dsub) {  
  ds <- dsub
  ## ds and dsub are the new variables; d is the original
  for (nm in names (ds)) {
    var <- sprintf ("d$%s", nm)
    A <- attributes (eval (parse (text=var)))
    if (!grepl ('Time', nm)) {
      A$dim[1] <- nrow(ds)
      A$class <- NULL
    } else {
      A$dim <- nrow (ds)
    }
    attributes (ds[,nm]) <- A
  }
  A <- attributes (d)
  A$Dimensions$Time$len <- nrow (ds)
  A$row.names <- 1:nrow (ds)
  A$names <- names (ds)
  attributes (ds) <- A
  return(ds)
}
for (Flight in 1:4) {
  standardFile <- sprintf ('%s%s/%srf%02dLAMS.nc', Directory, 
    Project, Project, Flight)
  MHFile <- sprintf ('%s%s/%srf%02d_LAMS_TAS_Update.nc', 
    Directory, Project, Project, Flight)
  netCDFfile = nc_open (standardFile)
  ATTG <- ncatt_get (netCDFfile, 0)   # get list of global attributes
  Time <- ncvar_get (netCDFfile, "Time")
  Time_units <- ncatt_get (netCDFfile, 'Time', 'units')
  ## special getNetCDFMatt.R uses this to get POSIX time
  Data <- getNetCDF(standardFile, 'ALL')
  DataMH <- getNetCDFMH(MHFile, 'ALL')
  if (Flight == 3) {Data <- Data[c(-1,-2),]}
  DataM <- cbind(Data, DataMH)
  # DataM <- transferAttributes (DataMH, DataM) # doesn't work
  DataM <- transferAttributes (Data, DataM)  ## from Data to DataM
  save(DataM, file=sprintf('DataM%d.Rdata', Flight))
}
## the above saves flight-1--4 data. 
## next: load and merge
load('DataM1.Rdata') 
DataM1 <- DataM 
DataM1$RF <- 1 
load('DataM2.Rdata') 
DataM2 <- DataM 
DataM2$RF <- 2 
load('DataM3.Rdata') 
DataM3 <- DataM 
DataM3$RF <- 3 
load('DataM4.Rdata') 
DataM4 <- DataM 
DataM4$RF <- 4 
N1 <- names(DataM1) 
N2 <- names(DataM2) 
N3 <- names(DataM3) 
N4 <- names(DataM4) 
DataM2 <- DataM2[, N2[N2 %in% N1]] 
DataM4 <- DataM4[, N4[N4 %in% N3]] 
## ensure that variable names are the same (need to exclude extraneous Time variable)
inx <- which(!(names(DataM4) %in% names(DataM1))) ## index to bad variable
DataM1 <- DataM1[, -inx] 
DataM2 <- DataM2[, -inx] 
DataM3 <- DataM3[, -inx] 
DataM4 <- DataM4[, -inx] 
DataM <- rbind(DataM1, DataM2) 
DataM <- rbind(DataM, DataM3) 
DataM <- rbind(DataM, DataM4) 
save(DataM, file='DataM.Rdata')
