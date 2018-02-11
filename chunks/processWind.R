
processWind <- function(Data) {
  Data$BEAM1speed <- zoo::na.approx (as.vector(Data$BEAM1speed), maxgap=1000, na.rm = FALSE)
  Data$BEAM2speed <- zoo::na.approx (as.vector(Data$BEAM2speed), maxgap=1000, na.rm = FALSE)
  Data$BEAM3speed <- zoo::na.approx (as.vector(Data$BEAM3speed), maxgap=1000, na.rm = FALSE)
  Data$BEAM4speed <- zoo::na.approx (as.vector(Data$BEAM4speed), maxgap=1000, na.rm = FALSE)
  ## replace with Matt's values
  if (PCA) {
    Data$BEAM1speed <- zoo::na.approx (as.vector(Data$V_LOS_Beam1), maxgap=1000, na.rm = FALSE)
    Data$BEAM2speed <- zoo::na.approx (as.vector(Data$V_LOS_Beam2), maxgap=1000, na.rm = FALSE)
    Data$BEAM3speed <- zoo::na.approx (as.vector(Data$V_LOS_Beam3), maxgap=1000, na.rm = FALSE)
    ## Beam4 is not processed by Matt; leave as SG value, but don't use below
    Data$BEAM4speed <- zoo::na.approx (as.vector(Data$BEAM4speed), maxgap=1000, na.rm = FALSE)
  }
  ## try to fix CTHDG_LAMS bad points in transition through 180: (may no longer be needed)
  if(FALSE) {
    for (i in 2:(nrow(Data)-1)) {
      if (is.na(Data$CTHDG_LAMS[i]) || is.na(Data$CTHDG_LAMS[i-1]) || is.na(Data$CTHDG_LAMS[i+1])) {next}
      if (abs(Data$CTHDG[i-1]-Data$CTHDG_LAMS[i+1]) > 10.) {next}
      if ((Data$CTHDG_LAMS[i-1] < 180.) && (Data$CTHDG_LAMS[i+1] > 180)) {Data$CTHDG_LAMS[i] <- NA}
      if ((Data$CTHDG_LAMS[i-1] > 180.) && (Data$CTHDG_LAMS[i+1] < 180)) {Data$CTHDG_LAMS[i] <- NA}
    }
  }
  Data$CTHDG_LAMS <- zoo::na.approx (as.vector(Data$CTHDG_LAMS), maxgap=1000, na.rm = FALSE)
  
  # also need the distances from the IRS to LAMS: (x,y,z)
  LL = c(-10.305, -6.319, 1.359)                # these are GV values
  # unit vectors along beams are then:
  #   a[i] = [cos(Theta[i]), -sin(Theta[i])*sin(Phi[i]), sin(Theta[i])*cos(Phi[i])]
  # and the dot products with the (i,j,k) unit vectors give the direction cosine matrix:
  S = c(cos(Theta[1]), -sin(Theta[1])*sin(Phi[1]), sin(Theta[1])*cos(Phi[1]), 
    cos(Theta[2]), -sin(Theta[2])*sin(Phi[2]), sin(Theta[2])*cos(Phi[2]), 
    cos(Theta[3]), -sin(Theta[3])*sin(Phi[3]), sin(Theta[3])*cos(Phi[3]))
  ## S4 is the four-beam solution, not used here
  S4 <- c(S, cos(Theta[4]), -sin(Theta[4])*sin(Phi[4]), sin(Theta[4])*cos(Phi[4]))
  dim(S) <- c(3,3)
  Si = t(solve(S))  # calculate the inverse of S -- this is the 3-beam version
  ## the following commented lines are the python code:
  # S4 = np.vstack ((S, [cos(Theta[3]), -sin(Theta[3])*sin(Phi[3]), sin(Theta[3])*cos(Phi[3])]))
  # StS =  linalg.inv (ma.dot (S4.T, S4))
  # M = ma.dot (StS, S4.T)      # matrix for finding relative wind from 4-beam LAMS
  dim(S4) <- c(3,4)
  StS <- S4 %*% t(S4)
  StS <- solve(StS)
  M <- StS %*% S4    ## this isn't used here
  
  A = c(Data$BEAM1speed, Data$BEAM2speed, Data$BEAM3speed)
  A4 <- c(A, Data$BEAM4speed)
  dim(A4) <- c(nrow(Data), 4)
  dim(A) <- c(nrow(Data), 3)
  RW = t (Si %*% t(A))    # gives u, v, w components, RW[,1] is u
  RW2 <- t (M %*% t(A4))
  ## calculate the error in the 4-beam solution:
  A4P <- t (t (S4) %*% t (RW2)) - A4
  CSQ <- A4P[,1]^2 + A4P[,2]^2 + A4P[,3]^2 + A4P[,4]^2
  
  ## set up a special data.frame for calculating wind:
  D <- data.frame("Time"=Data$Time)
  # D$TASX <- sqrt(RW2[,1]^2 + RW2[,2]^2 + RW2[,3]^2)  ## not used
  D$TASX3 <- sqrt(RW[,1]^2 + RW[,2]^2 + RW[,3]^2)    ## 3-beam version
  D$TASCSQ <- CSQ ## not used
  ## The following substitutions replace standard variables with LAMS-provided
  ## variables because the function WindProcessor() expects these names
  D$TASX <- D$TASX3   ## use the LAMS-provided value
  D$ATTACK <- atan (RW[, 3] / RW[, 1]) * 180 / pi
  D$SSLIP <-  atan (RW[, 2] / RW[, 1]) * 180 / pi
  D$GGVEW <- Data$CVEW_LAMS  
  D$GGVNS <- Data$CVNS_LAMS
  D$GGVSPD <- Data$CVSPD_LAMS
  D$VEW <- Data$CVEW_LAMS
  D$VNS <- Data$CVNS_LAMS
  D$THDG <- Data$CTHDG_LAMS
  D$ROLL <- Data$CROLL_LAMS
  D$PITCH <- Data$CPITCH_LAMS
  ## save some variables for transfer into DataM upon return:
  D$V_LOS_Beam1 <- Data$BEAM1speed
  D$V_LOS_Beam2 <- Data$BEAM2speed
  D$V_LOS_Beam3 <- Data$BEAM3speed
  D$V_LOS_Beam4 <- Data$BEAM4speed  ## this is really the SG-poly-provided value
  D$VXG <- RW[,1]
  D$VYG <- RW[,2]
  D$VZG <- RW[,3]
  # if (grepl('rf03', fname)) {
  #   rg <- setRange (D$Time, 192300,213000)
  #   D <- D[rg, ]
  # }
  ## use standard wind processor 
  ### rotation-rate corrections should be removed first, 
  ### but for standard routine they are insignificant so will be ignored for now.
  DW <- WindProcessor (data=D)
  return (DW)
}
