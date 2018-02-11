
SummarizeFit <- function(ft) {
  print (summary(ft)$call)
  print ("Coefficients:")
  print (summary(ft)$coefficients)
  print (sprintf ("Residual standard deviation: %.3f, dof=%d", summary(ft)$sigma, summary(ft)$df[2]))
  print (sprintf ("R-squared %.3f", summary(ft)$r.squared))
}
