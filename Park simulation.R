### Park Example ###
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

### synthetic function ###
park91a <- function(xx)
{
  if(xx[1]==0) xx[1] <- sqrt(.Machine$double.eps) # To prevent yielding infinity for Park function
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  
  term1a <- x1 / 2
  term1b <- sqrt(1 + (x2+x3^2)*x4/(x1^2)) - 1
  term1 <- term1a * term1b
  
  term2a <- x1 + 3*x4
  term2b <- exp(1 + sin(x3))
  term2 <- term2a * term2b
  
  y <- term1 + term2
  return(y)
}

park91alc <- function(xx)
{
  if(xx[1]==0) xx[1] <- sqrt(.Machine$double.eps) # To prevent yielding infinity for Park function
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  
  yh <- park91a(xx)
  
  term1 <- (1+sin(x1)/10) * yh
  term2 <- -2*x1 + x2^2 + x3^2
  
  y <- term1 + term2 + 0.5
  return(y)
}

### training data ###
n1 <- 40; n2 <- 20
d <- 4

rep <- 100
result.park.rmse <- matrix(NA, rep, 2)
result.park.meancrps <- matrix(NA, rep, 2)
result.park.comptime <- matrix(NA, rep, 2)
colnames(result.park.rmse) <- c("RNAmf", "Cokriging")
colnames(result.park.meancrps) <- c("RNAmf", "Cokriging") 
colnames(result.park.comptime) <- c("RNAmf", "Cokriging") 


for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  
  eps <- sqrt(.Machine$double.eps)
  
  ### test data ###
  x <- maximinLHS(1000, d)
  
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
                           lower=eps, #upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]

  
  result.park.rmse[i,1] <- sqrt(mean((predy-apply(x,1,park91a))^2)) 
  result.park.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,park91a))^2)) 
  
  result.park.meancrps[i,1] <- mean(crps(apply(x,1,park91a), predy, predsig2))
  result.park.meancrps[i,2] <- mean(crps(apply(x,1,park91a), pred.muficokm$mean, pred.muficokm$sig2))
  
  result.park.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.park.comptime[i,2] <- toc.cokm - tic.cokm
}

py_install("GPy")
py_install("pandas")
py_install("matplotlib")
py_install("scipy")
py_install("time")
py_run_file("python code/Park.py")
result.park.rmse <- cbind(result.park.rmse, NARGP=unlist(py$l2error))
result.park.meancrps <- cbind(result.park.meancrps, NARGP=unlist(py$meancrps))
result.park.comptime <- cbind(result.park.comptime, NARGP=unlist(py$comptime))
