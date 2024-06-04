### franke Example ###
franke2dl <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]

  term1 <- 0.75 * exp(-(9*x1-2)^2/4 - (9*x2-2)^2/4)
  term2 <- 0.75 * exp(-(9*x1+1)^2/49 - (9*x2+1)/10)
  term3 <- 0.5 * exp(-(9*x1-7)^2/4 - (9*x2-3)^2/4)
  term4 <- -0.2 * exp(-(9*x1-4)^2 - (9*x2-7)^2)

  y <- term1 + term2 + term3 + term4
  return(y)
}

franke2dm <- function(xx)
{
  y1 <- franke2dl(xx)
  y <- exp(-1.4*y1)*cos(3.5*pi*y1)
  return(y)
}

franke2dh <- function(xx)
{
  y2 <- franke2dm(xx)
  y <- sin(2*pi*(y2-1))
  return(y)
}


### training data ###
n1 <- 20; n2 <- 15; n3 <- 10; rep <- 100
d <- 2
### test data ###
n.test <- 1000

cat("Franke function: n1 =", n1, ", n2 =", n2, ", n3 =", n3, ", rep =", rep, "\n")


result.franke.rmse <- matrix(NA, rep, 3)
result.franke.meancrps <- matrix(NA, rep, 3)
result.franke.comptime <- matrix(NA, rep, 3)
colnames(result.franke.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.franke.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.franke.comptime) <- c("RNAmf", "Cokriging", "NARGP")

for(i in 1:rep) {
  set.seed(i)

  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  X3 <- maximinLHS(n3, d)

  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))

  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)

  y1 <- apply(X1,1,franke2dl)
  y2 <- apply(X2,1,franke2dm)
  y3 <- apply(X3,1,franke2dh)

  ### test data ###
  x <- maximinLHS(n.test, d)
  
  # Generate the filename using paste0
  filename <- paste0("RDSfile/file", i, ".rds") # change path
  saveRDS(list(X1=X1, X2=X2, X3=X3, Y1=y1, Y2=y2, Y3=y3, Xtest=x, Ytest=apply(x,1,franke2dh)), file = filename)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy <- pred.RNAmf$mu
  predsig2 <- pred.RNAmf$sig2
  toc.RNAmf <- proc.time()[3]

  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.5,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]


  result.franke.rmse[i,1] <- sqrt(mean((predy-apply(x,1,franke2dh))^2))
  result.franke.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,franke2dh))^2))
  # result.franke.rmse[i,3] <- py$error

  result.franke.meancrps[i,1] <- mean(crps(apply(x,1,franke2dh), predy, predsig2))
  result.franke.meancrps[i,2] <- mean(crps(apply(x,1,franke2dh), pred.muficokm$mean, pred.muficokm$sig2))
  # result.franke.meancrps[i,3] <- py$crps

  result.franke.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.franke.comptime[i,2] <- toc.cokm - tic.cokm
  # result.franke.comptime[i,3] <- py$ctime

  # boxplot(result.franke.rmse[1:i,, drop=FALSE])
}

