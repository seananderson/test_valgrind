#' @title fit linreg
#'
#' @description fits linreg.
#'
#' @importFrom TMB MakeADFun sdreport
#' @importFrom mgcv gam
#'
#' @param x independent.
#' @param y dependent.
#'
#' @examples
#' fit()
#'
#' @useDynLib testsmoothers, .registration = TRUE
#' @export
fit <-
  function(){

    # library(TMB)
    # library(mgcv) #Use gam
    library(Matrix) #Use sparse matrices

    # Tutorial: sparse matrixes in R and TMB
    M = matrix(1:4,2,2)   # ordinary 2x2 matrix
    M_block_diag = .bdiag(list(M,M)) # Block diagonal (sparse) matrix
    # data.class(M_block_diag)     # Check data.class
    # print(M_block_diag)    # dots means 0 value

    #Load the data and compile c++ code------
    # Vegetation <- read.table(file = "Vegetation.txt", header = TRUE, dec = ".")
    # Vegetation = Vegetation[!is.na(Vegetation$Richness),]
    Vegetation <- structure(list(SAMPLEYR = c(1958L, 1962L, 1967L, 1974L, 1981L,
          1994L, 2002L, 1958L, 1962L, 1967L, 1974L, 1981L, 1994L, 2002L,
          1958L, 1962L, 1967L, 1974L, 1981L, 1989L, 1994L, 2002L, 1958L,
          1962L, 1967L, 1974L, 1981L, 1989L, 1994L, 2002L, 1958L, 1962L,
          1967L, 1974L, 1981L, 1989L, 1994L, 2002L, 1958L, 1962L, 1967L,
          1974L, 1981L, 1989L, 1994L, 2002L, 1962L, 1967L, 1974L, 1981L,
          1994L, 2002L, 1962L, 1967L, 1974L, 1981L, 1994L, 2002L), Time = c(1L,
          2L, 3L, 4L, 5L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 7L, 8L, 1L, 2L, 3L,
          4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L,
          4L, 5L, 6L, 7L, 8L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 2L, 3L, 4L,
          5L, 7L, 8L, 2L, 3L, 4L, 5L, 7L, 8L), Transect = c(1L, 1L, 1L,
          1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L,
          3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L,
          5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L,
          7L, 8L, 8L, 8L, 8L, 8L, 8L), Richness = c(8L, 6L, 8L, 8L, 10L,
          7L, 6L, 5L, 8L, 6L, 6L, 6L, 6L, 6L, 7L, 10L, 8L, 18L, 12L, 11L,
          7L, 10L, 8L, 9L, 6L, 12L, 13L, 10L, 8L, 8L, 13L, 16L, 9L, 14L,
          11L, 13L, 11L, 12L, 9L, 10L, 14L, 14L, 10L, 14L, 9L, 12L, 11L,
          12L, 14L, 9L, 5L, 12L, 9L, 10L, 16L, 12L, 10L, 14L), ROCK = c(27,
          26, 30, 18, 23, 26, 39, 25, 24, 21, 18, 19, 10.5, 26, 56, 45,
          28, 27, 10, 17, 26.5, 36, 53, 59, 45, 41, 35, 44, 53.5, 59, 15,
          20, 4, 8, 5, 4, 2, 8, 5, 7, 5, 6, 3, 2, 3, 0, 20, 7, 11, 6, 0,
          15, 25, 23, 14, 11, 8, 13), LITTER = c(30, 20, 24, 35, 22, 26,
          19, 26, 24, 16, 25, 28, 41.5, 18, 17, 7, 14, 15, 37, 17, 14,
          19, 10, 5, 9, 12, 24, 10, 18, 9, 23, 21, 51, 34, 28, 30, 32,
          29, 32, 20, 29, 19, 23, 32, 22.5, 28, 26, 29, 23, 40, 14.5, 21,
          24, 15, 15, 18, 29, 26), BARESOIL = c(26, 28, 30, 16, 9, 23,
          19, 33, 29, 41, 33, 14, 29, 42, 16, 23, 37, 7, 0, 33, 17.5, 4,
          20, 13, 30, 7, 3, 12, 7.5, 5, 20, 19, 13, 2, 1, 11, 11, 5, 30,
          37, 12, 13, 9, 4, 16.5, 17, 27, 26, 0, 0, 35, 7, 20, 20, 14,
          4, 27, 13), FallPrec = c(30.21999931, 99.55999756, 43.43000031,
          54.86000061, 24.37999916, 10.15999985, 34.29000092, 30.21999931,
          99.55999756, 43.43000031, 54.86000061, 24.37999916, 10.15999985,
          34.29000092, 33.77999878, 143.5099945, 56.38000107, 68.31999969,
          57.40000153, 17.52000046, 21.07999992, 29.20999908, 33.77999878,
          143.5099945, 56.38000107, 68.31999969, 57.40000153, 17.52000046,
          21.07999992, 29.20999908, 41.13999939, 153.4100037, 75.43000031,
          72.63999939, 41.65000153, 9.899999619, 32, 48.50999832, 41.13999939,
          153.4100037, 75.43000031, 72.63999939, 41.65000153, 9.899999619,
          32, 48.50999832, 153.4100037, 75.43000031, 72.63999939, 41.65000153,
          32, 48.50999832, 153.4100037, 75.43000031, 72.63999939, 41.65000153,
          32, 48.50999832), SprTmax = c(15.77000046, 15.21000004, 12.76000023,
          14, 14.32999992, 16.90999985, 13.85999966, 15.78999996, 15.22999954,
          12.77999973, 14.01000023, 14.35000038, 16.92000008, 13.88000011,
          11.78999996, 11.25, 8.909999847, 10.09000015, 10.42000008, 11,
          12.96000004, 9.880000114, 11.92000008, 11.38000011, 9.039999962,
          10.22000027, 10.55000019, 11.14000034, 13.09000015, 10.01000023,
          11.86999989, 11.43999958, 9.25, 10.35999966, 10.68999958, 11.22000027,
          13.15999985, 9.899999619, 12.02999973, 11.59000015, 9.399999619,
          10.51000023, 10.84000015, 11.36999989, 13.31000042, 10.05000019,
          12.81999969, 10.56000042, 11.69999981, 12.02999973, 14.56000042,
          11.35000038, 12.93000031, 10.67000008, 11.81000042, 12.14999962,
          14.67000008, 11.46000004)), row.names = c(1L, 2L, 3L, 4L, 5L,
        7L, 8L, 9L, 10L, 11L, 12L, 13L, 15L, 16L, 17L, 18L, 19L, 20L,
        21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 33L,
        34L, 35L, 36L, 37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L,
        47L, 48L, 50L, 51L, 52L, 53L, 55L, 56L, 58L, 59L, 60L, 61L, 63L,
        64L), class = "data.frame")

    #Set up spline structure by using mgcv---
    gam_setup = gam(Richness ~ s(ROCK, bs = "cs") +
      s(LITTER, bs = "cs") + s(BARESOIL, bs = "cs") +
      s(FallPrec, bs = "cs") + s(SprTmax, bs = "cs"),
    data = Vegetation,fit=FALSE)

    #Extrtact penelization matrices
    S_ROCK = gam_setup$smooth[[1]]$S[[1]]
    S_LITTER = gam_setup$smooth[[2]]$S[[1]]
    S_BARESOIL = gam_setup$smooth[[3]]$S[[1]]
    S_FallPrec = gam_setup$smooth[[4]]$S[[1]]
    S_SprTmax = gam_setup$smooth[[5]]$S[[1]]

    S_list = list(S_ROCK,S_LITTER,S_BARESOIL,S_FallPrec,S_SprTmax)
    S_combined = .bdiag(S_list)         # join S's in sparse matrix
    Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S
    #----------------------------------------

    #For report, used for constructing plots----
    ROCK=seq(min(Vegetation$ROCK),max(Vegetation$ROCK),by = 1)
    LITTER=seq(min(Vegetation$LITTER),max(Vegetation$LITTER),by = 1)
    BARESOIL=seq(min(Vegetation$BARESOIL),max(Vegetation$BARESOIL),by = 1)
    FallPrec=seq(min(Vegetation$FallPrec),max(Vegetation$FallPrec),by = 0.2)
    SprTmax=seq(min(Vegetation$SprTmax),max(Vegetation$SprTmax),by = 0.2)

    rockReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(ROCK))
    litterReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(LITTER))
    soilReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(BARESOIL))
    fallReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(FallPrec))
    sprReport = PredictMat(gam_setup$smooth[[5]],data = data.frame(SprTmax))

    designMatrixForReport = list(rockReport,litterReport,soilReport,fallReport,sprReport)
    #-------------------------------------------

    #Define data object which is given to TMB---
    data = list(Y = Vegetation$Richness, # Response
      X = gam_setup$X[,-1],  # Design matrix, without intercept
      S = S_combined,      # Combined penalty matrix
      Sdims = Sdims,
      designMatrixForReport = .bdiag(designMatrixForReport),
      flag = 1)
    #-------------------------------------------

    #Define parameter object given to TMB-------
    par = list(
      beta0 = 0,  # Intercept
      beta = rep(0,sum(Sdims)),  # Spline coefficients
      log_lambda = rep(rep(0,length(Sdims))), #Log spline penalization coefficients
      log_sigma = 0
    )
    #-------------------------------------------

    #Fit model----------------------------------
    obj = MakeADFun(data = data, parameters = par,random="beta",DLL = "testsmoothers", silent = TRUE)
    #obj <- normalize(obj, flag="flag")
    opt = nlminb(obj$par,obj$fn,obj$gr)
    rep = sdreport(obj)
    #-------------------------------------------

    rep
  }
