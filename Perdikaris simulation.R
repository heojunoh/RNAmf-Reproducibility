### Perdikaris Example ###
### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}

### training data ###
n1 <- 13; n2 <- 8

rep <- 100
result.perd.rmse <- matrix(NA, rep, 2)
result.perd.meancrps <- matrix(NA, rep, 2)
result.perd.comptime <- matrix(NA, rep, 2)
colnames(result.perd.rmse) <- c("RNAmf", "Cokriging")
colnames(result.perd.meancrps) <- c("RNAmf", "Cokriging") 
colnames(result.perd.comptime) <- c("RNAmf", "Cokriging") 

for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  
  eps <- sqrt(.Machine$double.eps)
  
  ### test data ###
  x <- seq(0,1,length.out=1000)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predRNAmf(fit.RNAmf, x)
  predy <- pred.RNAmf$mu
  predsig2 <- pred.RNAmf$sig2
  toc.RNAmf <- proc.time()[3]

  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  
  result.perd.rmse[i,1] <- sqrt(mean((predy-f2(x))^2)) 
  result.perd.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-f2(x))^2)) 
  
  result.perd.meancrps[i,1] <- mean(crps(f2(x), predy, predsig2))
  result.perd.meancrps[i,2] <- mean(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) 
  
  result.perd.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.perd.comptime[i,2] <- toc.cokm - tic.cokm
}

py_run_file("python code/Perdikaris.py")
result.perd.rmse <- cbind(result.perd.rmse, NARGP=unlist(py$l2error))
result.perd.meancrps <- cbind(result.perd.meancrps, NARGP=unlist(py$meancrps))
result.perd.comptime <- cbind(result.perd.comptime, NARGP=unlist(py$comptime))
