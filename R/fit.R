#' @title fit linreg
#'
#' @description fits linreg.
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom stats optim
#'
#' @param x independent.
#' @param y dependent.
#'
#' @examples
#' x = 1:10
#' y = rnorm(10) + x
#' out = fit( x=x, y=y )
#'
#' @useDynLib testvalgrind, .registration = TRUE
#' @export
fit <-
function( x, y ){
  data <- list(Y = y, x=x)
  parameters <- list(a=0, b=0, logSigma=0)
  obj <- MakeADFun(data, parameters, DLL="testvalgrind")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  #opt
  #opt$hessian ## <-- FD hessian from optim
  #obj$he()    ## <-- Analytical hessian
  sdrep = sdreport(obj)

  out = list(
    obj = obj,
    opt = opt,
    sdrep = sdrep
  )
  return(out)
}
