## ----initialization, include=FALSE---------------------------------------

require(knitr)
opts_chunk$set(fig.path='figure/SO-', echo=FALSE, include=FALSE, fig.lp="fig:", dev='png', dpi=100, fig.show='hold', size='footnotesize', replace.assign=TRUE, width=49)
opts_chunk$set(fig.width=6, fig.height=5, fig.align="center", digits=4)
options(digits=5)
thisFileName <- "WindInSOCRATES"
require(Ranadu, quietly = TRUE, warn.conflicts=FALSE)
require(ggplot2)
require(grid)
require(ggthemes)
require(zoo)
library(scales)
source('chunks/VSpec.R')  ## temporary, pending inclusion in Ranadu
source('chunks/removeSpikes.R')
source('chunks/DemingFit.R')    ## temporary, pending Ranadu update
source('chunks/SplitDV.R')
refline <- function (vmin=-100, vmax=100) {
  lines(c(vmin, vmax), c(vmin, vmax), col='darkorange', lwd=2, lty=2)
}
ReviseProjects <- c('SOCRATES')  ## these are the projects to process 

Directory <- DataDirectory()

## ----nVarCalc, eval=FALSE------------------------------------------------
## 
Directory <- DataDirectory ()
source('chunks/AddWind.R')
for (Project in ReviseProjects) {
  ## get the list of flights, low rate:
  Fl <- sort (list.files (sprintf ("%s%s/", Directory, Project),
    sprintf ("%srf...nc$", Project)))
  for (fn in Fl[1]) {
    fname <- sprintf('%s%s/%s', Directory, Project, fn)
    fnew <- sub ('.nc', 'Y.nc', fname)
    Z <- file.copy (fname, fnew, overwrite=TRUE)  ## BEWARE: overwrites without warning!!
    ## read variables needed for the calculation
    FI <- DataFileInfo (fname, LLrange=FALSE)
    
    ## for some old projects:
    if (!('GGVSPD' %in% FI$Variables)) {
      if ('GGVSPDB' %in% FI$Variables) {
        VR [which (VR == 'GGVSPD')] <- 'GGVSPDB'
      } else if ('VSPD_A' %in% FI$Variables) {
        VR [which (VR == 'GGVSPD')] <- 'VSPD_A'
      } else if ('VSPD_G' %in% FI$Variables) {
        VR [which (VR == 'GGVSPD')] <- 'VSPD_G'
      } else {
        print ('ERROR: no VSPD variable found')
        exit()
      }
    }
    for (Var in VR) {
      if (!(Var %in% FI$Variables)) {
        print (sprintf (' required variable %s not found in file %s; skipping...', Var, fname))
        exit()
      }
    }
    ## 
    DY <- getNetCDF(fname, VR)
   
    
    ## ----newNetCDF, eval=FALSE-----------------------------------------------
    ## 
    source ('chunks/copyAttributes.R')
    ## 
    netCDFfile <- nc_open (fnew, write=TRUE)
    Dimensions <- attr (Data, "Dimensions")
    Dim <- Dimensions[["Time"]]
    Rate <- 1
    if ("sps25" %in% names (Dimensions)) {
      Rate <- 25
      Dim <- list(Dimensions[["sps25"]], Dimensions[["Time"]])
    }
    if ("sps50" %in% names (Dimensions)) {
      Rate <- 50
      Dim <- list(Dimensions[["sps50"]], Dimensions[["Time"]])
    }
    addAKY <- TRUE
    addGP <- FALSE
    addTC <- TRUE
    addROC <- TRUE
    Data <- AddWind(DY, Rate, addAKY, addGP, addTC, addROC)    ## default adds everything
    ## 
    DATT <- Data  ## save to ensure that attributes are preserved
    
    ## variables to add to the netCDF file:
    VarNew <- c('AKY', 'WIY', 'QCTC', 'AK_GP', 'SS_GP', 'WIG', 'WDG', 'WSG', 'TASG', 'UXG', 'VYG', 'ROC', 'TASTC', 'WDTC', 'WSTC', 'WITC', 'UXTC', 'VYTC')
    VarOld <- c('AKRD', 'WIC', 'QCFC', 'AKRD', 'SSRD', 'WIC', 'WDC', 'WSC', 'TASX', 'UXC', 'VYC', 'GGVSPD', 'TASX', 'WDC', 'WSC', 'WIC', 'UXC', 'VYC')
    VarUnits <- c('degrees', 'm/s', 'degrees', 'degrees', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'degrees', 'm/s', 'm/s', 'm/s', 'm/s')
    VarStdName <- c('angle-of-attack, CF', 'vertical wind, CF', 'dynamic pressure, pitot-static, corrected',
      'angle-of-attack, GP', 'sideslip angle, GP', 
      'vertical wind, GP', 'wind direction, GP', 'wind speed, GP', 'true airspeed, GP', 
      'wind longitudinal component, GP', 'wind lateral component, GP', 'rate of climb', 
      'true airspeed, pitot-static', 'wind direction, pitot-static', 'wind speed, pitot-static', 
      'vertical wind, pitot-static', 'wind longitudinal component, pitot-static', 'wind lateral component, pitot-static')
    VarLongName <- c('angle of attack, complementary-filter',
      'vertical wind using comp-filter angle of attack',
      'dynamic pressure from the pitot-static sensor, corrected',
      'angle of attack from the gustpod',
      'sideslip angle from the gustpod',
      'vertical wind from the gustpod',
      'horizontal wind direction from the gustpod',
      'horizontal wind speed from the gustpod',
      'true airspeed from the gustpod',
      'horizontal wind, longitudinal component, gustpod',
      'horizontal wind, lateral component, gustpod',
      'rate of climb of the aircraft from pressure',
      'true airspeed from the pitot-static sensor',
      'wind direction based on the pitot-static airspeed',
      'wind speed based on the pitot-static airspeed',
      'vertical wind based on TASTC and AKY',
      'horizontal wind, longitudinal component, pitot-static',
      'horizontal wind, lateral component, pitot-static')
    
    ## create the new variables
    varCDF <- list ()
    for (i in 1:length(VarNew)) {
      if (!addAKY && (i <= 2)) {next}
      if (!addGP && (i %in% 4:11)) {next}
      if (!addTC && (i %in% c(3, 13:18))) {next}
      if (!addROC && (i == 12)) {next}
      print (sprintf ('new-netcdf %d%% done', as.integer(100*(i-1)/length(VarNew))))
      varCDF[[i]] <- ncvar_def (VarNew[i],
        units=VarUnits[i],
        dim=Dim,
        missval=as.single(-32767.), prec='float',
        longname=VarLongName[i])
      if (i == 1) {
        newfile <- ncvar_add (netCDFfile, varCDF[[i]])
      } else {
        newfile <- ncvar_add (newfile, varCDF[[i]])
      }
      ATV <- ncatt_get (netCDFfile, VarOld[i])
      copy_attributes (ATV, VarNew[i], newfile)
      ncatt_put (newfile, VarNew[i], attname="standard_name",
        attval=VarStdName[i])
      if (Rate == 1) {
        ncvar_put (newfile, varCDF[[i]], Data[, VarNew[i]])
      } else if (Rate == 25) {
        ncvar_put (newfile, varCDF[[i]], Data[, VarNew[i]], count=c(25, nrow(Data)/25))
      }
    }
    nc_close (newfile)
  }
}


## ----make-zip-archive, INCLUDE=TRUE, eval=FALSE--------------------------
## 
## cat (toLatex(sessionInfo()), file="SessionInfo")
## system (sprintf("zip WindInSOCRATES.zip WindInSOCRATES.Rnw WindInSOCRATES.pdf NoteReSOCRATESwindProcessing.pdf WorkflowWindInSocrates.pdf WAC.bib ./chunks/* SessionInfo"))
## 

