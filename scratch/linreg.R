#library(TMB)
#compile("linreg.cpp")
#dyn.load(dynlib("linreg"))

library(testvalgrind)
set.seed(123)
x = 1:10
y = rnorm(10) + x

#fit <- function( x, y ){
#  data <- list(Y = y, x=x)
#  parameters <- list(a=0, b=0, logSigma=0)
#  obj <- MakeADFun(data, parameters, DLL="linreg")
#  obj$hessian <- TRUE
#  opt <- do.call("optim", obj)
#  #opt
#  #opt$hessian ## <-- FD hessian from optim
#  #obj$he()    ## <-- Analytical hessian
#  sdrep = sdreport(obj)
#
#  out = list(
#    obj = obj,
#    opt = opt,
#    sdrep = sdrep
#  )
#  return(out)
#}
#

fit( x=x, y=y )
