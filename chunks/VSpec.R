#' @title VSpec
#' @description Produces a variance spectrum in a ggplot2 file suitable for printing.
#' @details For the variable provided, which must be in the supplied data.frame 
#' that must also contain the variables "Time" and "TASX' and a "Rate" attribute, this 
#' function constructs a plot of the spectral variance of that variable using either 
#' the standard "spectrum" function of R (the default) or the maximum-entropy
#' functions of Ranadu. Some options are provided for each, with defaults that are usually
#' appropriate. If the ADD parameter is set to the ggplot definition returned from a
#' previous call, the new spectrum is added to that plot (up to four plotted spectra
#' on a single plot). The mean and trend are removed before constructing the spectrum.
#' Optional additional smoothing in logarithmic intervals in frequency is available
#' through the "smoothBins" parameter. The plot includes background lines showing the
#' -5/3 slope expected for homogeneous isotropic turbulence, and for wind variables
#' the magnitude of these lines represent factor-of-10 changes in the eddy dissipation
#' rate with the larger-dot line representing 1e-4 m^2/s^3. 
#' @aliases vSpec vspec
#' @author William Cooper
#' @import scales
#' @importFrom zoo na.approx
#' @export VSpec
#' @param .data A data.frame containing at least the variables "Time" and ".Variable" where
#' ".Variable" is the second (required) parameter. It should also have an attribute "Rate"
#' if its rate is different from 1 Hz (the default). Any restrictions on the time range
#' should be applied to the data.frame before it is supplied to this function. See the
#' examples below. If subsetting removes the "Rate" attributes from the data.frame, the
#' value from the original data.frame should be added to .data.
#' @param .Variable The name of a variable that is a column in .data and for which the
#' variance spectrum will be constructed. (Otherwise an error message is generated.)
#' @param VLabel A character string that will be used as the label for this variable in
#' the legend. The default is .Variable.
#' @param col The color to use when plotting this variable. The default is NA, and
#' in this case the following plot colors will be used in order: blue, forestgreen, red, 
#' azure.
#' @param type Two choices are avaiable at present: "spectrum" (the default), which
#' uses the R routine "spectrum()" from the stats package to estimate the spectral
#' density, and "MEM" (or any value other than "spectrum") to use the maximum-entropy
#' method of spectral estimation as implemented in the Ranadu routines memCoef() and
#' memEstimate(). The parameter "spans" applies only to the "spectrum" method, and the
#' arguments "poles" and "resolution" apply only to the "MEM" method.
#' @param method This is the same as "type" and over-rides it if present.
#' @param spans An odd integer (or forced odd if even) specifying the number of frequencies
#' to span when averaging the spectral variance estimate produced by the R routine
#' "spectrom". See help for that function for more information about the nature of
#' this averaging. The default value is 49. If spans=NULL or spans <= 4 this averaging
#' is suppressed.
#' @param ae A factor that scales the predicted lines representing constant eddy dissipation
#' rate to adjust for longitudinal (where ae should be 0.15) or lateral (0.2) spectra. The
#' default is 0.2; use 0.15 for the variables TASX, UXC, etc. This applies to only the first
#' spectrum plotted; others use the same eddy-dissipation-rate scaling. To make an
#' appropriate scale adjustment for plots that include both longitudinal and lateral
#' variance spectra, consider multiplying the longitudinal variable by sqrt(0.15/0.2).
#' @param smoothBins If a value larger than 5 is provided, the frequency range is divided 
#' into this number of intervals evenly spaced in the logarithm of the frequency. Then
#' estimates of the spectral density are binned into those intervals and averaged to
#' smooth the spectrum. Initial smoothing can be provided by "spans" (if larger than 4) for
#' the "spectrum" method and by using a small number of poles for the "MEM" method; the
#' smoothing by the "smoothBins" parameter is applied after and in addition to those
#' smoothing methods. The default (0) suppresses this smoothing.
#' @param poles The number of poles to use for the maximum-entropy estimates. The default (50)
#' usually is a good first choice. For more structure, try 100. If you go beyond around 150,
#' the method becomes too slow to be practical, esp. for high-rate files.
#' @param resolution The increment (as a fraction of the logarithm of the frequency range) 
#' at which to evaluate the MEM estimate of spectral variance If a fine feature is expected
#' in the spectrum, the resolution should be small enough to isolate it, but small values
#' increase the variance in the individual estimates.
#' @param ADD This parameter has the default value NA, which causes the function to plots 
#' only the spectrum for the variable provided. If a spectrum for an other set of variables 
#' has already been defined by previous calls to VSpec, setting ADD to the plot definition 
#' returned by that previous call will add this plot to the previous plot. Up to four
#' variables can be included in the final plot. See the examples.
#' @return A ggplot2 definition for the plot of spectral density as a function of frequency.
#' The normalization is one-sided; i.e., the integral of the spectral variance from zero
#' to infinity is the total variance of the variable. The resulting plot definition
#' can be plotted (via, e.g., 'print (VarSpec(...))) or
#' saved for later addition of more variables or for later plotting. The plot is returned
#' with the standard ggplot theme; to use the Ranadu theme "theme_WAC()", add it to the
#' plot definition that is returned before plotting. In addition, to make it possible to 
#' superimpose future plots, the following variables are saved in the global environment: 
#' .clinesVSpec and .VSpecDF{1,2,3}, .VSpecVar{2,3}. The function does not check for 
#' collision with other possible uses of those names in the global environment.
VSpec <- function (.data, .Variable, VLabel=NA, col=NA, type='spectrum', method=NA, 
                   spans=49, ae=0.2, smoothBins=0, poles=50, resolution=0.0001, ADD=NA) {
  if (is.data.frame(.data)) {
    if (.Variable %in% names(.data)) {
      Z <- capture.output (v <- detrend (.data[, c('Time', .Variable)]))
      if (!is.na(VLabel)) {    ## use this alternate name in legend
        V <- VLabel
      } else {
        V <- .Variable
      }
    } else {
      print(sprintf('VSpec ERROR: Variable %s is not in the supplied data.frame', .Variable))
      return (NA)
    }
  } else {
    print('VSpec ERROR: first argument is not a data.frame.')
    return (NA)
  }
  if (is.null(attr(.data, 'Rate'))) {
    print ('VSpec warning: Rate attribute missing from data.frame, so using Rate=1')
    Rate <- 1
  } else {
    Rate <- attr(.data, 'Rate')
  }
  if (!is.na(method)) {type <- method}  ## method over-rides if present
  if (type == 'spectrum') {
    if (!is.null(spans)) {
      if (!(spans %% 2)) {spans <- spans + 1}
      if (spans <= 5) {spans <- NULL}
    }
    S <- spectrum (ts(SmoothInterp(v, .Length=0), frequency=Rate), span=spans, plot=FALSE)
    freq <- S$freq
    fpf <- 2 * S$spec * freq
  } else {  ## MEM section
    A <- memCoef (v, poles)
    ld <- nrow(.data)
    fmin <- log (Rate / ld)
    fmax <- log (0.5*Rate)
    bins <- as.integer (1/resolution)
    df <- (fmax-fmin) / bins
    fdtl <- fmin + df * (0:bins)
    freq <- exp (fdtl)
    psComplex <- memEstimate (freq / Rate, A) / Rate
    ps <- 2 * Rate * Mod (psComplex)^2
    fpf <- freq * ps
  }
  
  if(smoothBins > 9) {
    bs1 <- binStats(data.frame(fpf, log(freq)), bins=smoothBins)
    bs1 <- rbind (bs1, data.frame(xc=bs1$xc[nrow(bs1)], ybar=bs1$ybar[nrow(bs1)],
      sigma=bs1$sigma[nrow(bs1)], nb=1))
    freq <- exp(bs1$xc)
    fpf <- bs1$ybar
    bs1$sigma <- ifelse (bs1$nb > 2, bs1$sigma/sqrt(bs1$nb), NA)
    rna <- is.na(bs1$sigma)
    bs1$sigma[rna] <- bs1$ybar[rna] / 2
    # bs1 <<- bs1
  }
  
  DF <- data.frame(freq, fpf)
  if (is.na(ADD[1])) {
    ## first call: redefine VSpecDF
    assign('.VSpecDF1', DF, envir=.GlobalEnv)
    labx <- 'frequency [Hz]'
    laby <- 'fP(f)'
    xlim <- c(0.001,15)
    ylim <- c(0.0001, 1) 
    g <- ggplot(data=.VSpecDF1)         
    g <- g + geom_path (aes(x=freq, y=fpf, colour=V), na.rm=TRUE) +  
             xlab(labx) + ylab (laby) 
    if (is.na(col)) {col <- 'blue'}
    .clinesVSpec <- col
    names(.clinesVSpec) <- V
    .clinesVSpec <<- .clinesVSpec    
  } else {
    ## assign name based on elements in clinesVSpec
    N <- length(.clinesVSpec) + 1
    if (is.na(col)) {
      if (N == 2) {col <- 'forestgreen'}
      if (N == 3) {col <- 'brown'}
      if (N == 4) {col <- 'cyan'}
    }
    nc <- names(.clinesVSpec)
    .clinesVSpec <- c(.clinesVSpec, col)
    names(.clinesVSpec) <- c(nc, V)
    .clinesVSpec <<- .clinesVSpec
    VName <- sprintf('.VSpecDF%d', N)
    assign(VName, DF, pos=.GlobalEnv)
    if (N == 2) {
      .VSpecVar2 <<- V
      g <- ADD + geom_path (aes(x=freq, y=fpf, colour=.VSpecVar2), data=get(VName), na.rm=TRUE)
    } else if (N == 3) {
      .VSpecVar3 <<- V
      g <- ADD + geom_path (aes(x=freq, y=fpf, colour=.VSpecVar3), data=get(VName), na.rm=TRUE)
    } else if (N == 4) {
      .VSpecVar4 <<- V
      g <- ADD + geom_path (aes(x=freq, y=fpf, colour=.VSpecVar4), data=get(VName), na.rm=TRUE)  
    }
  }
  g <- suppressMessages(g + scale_colour_manual (name='', values=.clinesVSpec))
  if (is.na(ADD[1])) {
    g <- g + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n=4), #limits = xlim, 
      labels = trans_format("log10", math_format(10^.x))) +           
      scale_y_log10(breaks =            trans_breaks("log10", function(x) 10^x, n=4), #limits = ylim,             
        labels = trans_format("log10", math_format(10^.x))) +           
      coord_cartesian(xlim=xlim, ylim=ylim) 
    tasAverage <- mean(.data$TASX, na.rm=TRUE)
    for (i in (-8:0)) {
      a = ae * 10.^(i*(2/3)) * tasAverage^(2/3)
      lw = ifelse(i == -4, 1.2, 0.5)
      DFL <- data.frame(x=xlim, y=c(a/xlim[1]^(2/3), a/xlim[2]^(2/3)))
    # print(DFL)
      g <- g + geom_path (data=DFL, aes(x=x, y=y), colour='darkorange', lwd=lw, lty=3)
    }
    # g <- g + theme_WAC()
  }
  
  return(g)
}

