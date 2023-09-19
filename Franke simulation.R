### Franke Example ###
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

### synthetic function ###
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
n1 <- 20; n2 <- 15; n3 <- 10

rep <- 100
result.Franke.rmse <- matrix(NA, rep, 2)
result.Franke.meancrps <- matrix(NA, rep, 2)
result.Franke.comptime <- matrix(NA, rep, 2)
colnames(result.Franke.rmse) <- c("RNAmf", "Cokriging")
colnames(result.Franke.meancrps) <- c("RNAmf", "Cokriging")
colnames(result.Franke.comptime) <- c("RNAmf", "Cokriging")

for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  X1 <- maximinLHS(n1, 2)
  X2 <- maximinLHS(n2, 2)
  X3 <- maximinLHS(n3, 2)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  y1 <- apply(X1,1,franke2dl)
  y2 <- apply(X2,1,franke2dm)
  y3 <- apply(X3,1,franke2dh)
  
  
  eps <- sqrt(.Machine$double.eps)

  ### test data ###
  x <- maximinLHS(1000, 2)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predRNAmf2(fit.RNAmf, x)
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
  
  
  result.Franke.rmse[i,1] <- sqrt(mean((predy-apply(x,1,franke2dh))^2)) 
  result.Franke.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,franke2dh))^2)) 
  
  result.Franke.meancrps[i,1] <- mean(crps(apply(x,1,franke2dh), predy, predsig2)) 
  result.Franke.meancrps[i,2] <- mean(crps(apply(x,1,franke2dh), pred.muficokm$mean, pred.muficokm$sig2)) 
  
  result.Franke.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.Franke.comptime[i,2] <- toc.cokm - tic.cokm
}

py_install("GPy")
py_install("pandas")
py_install("matplotlib")
py_install("scipy")
py_install("time")
py_run_file("python code/Franke.py")
result.Franke.rmse <- cbind(result.Franke.rmse, NARGP=unlist(py$l2error))
result.Franke.meancrps <- cbind(result.Franke.meancrps, NARGP=unlist(py$meancrps))
result.Franke.comptime <- cbind(result.Franke.comptime, NARGP=unlist(py$comptime))
