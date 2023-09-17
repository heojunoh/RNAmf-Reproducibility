### Perdikaris Example ###
library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(RNAmf)
library(reticulate)

crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

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
result.perd.rmse2 <- matrix(NA, rep, 2)
result.perd.meancrps <- matrix(NA, rep, 2)
result.perd.comptime <- matrix(NA, rep, 2)
colnames(result.perd.rmse) <- c("closed", "Cokriging")
colnames(result.perd.meancrps) <- c("closed", "Cokriging") 
colnames(result.perd.comptime) <- c("closed", "Cokriging") 

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
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.closed <- predRNAmf(fit.closed, x)
  predy <- pred.closed$mu
  predsig2 <- pred.closed$sig2
  toc.closed <- proc.time()[3]

  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  
  result.perd.rmse[i,1] <- sqrt(mean((predy-f2(x))^2)) # closed form
  result.perd.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-f2(x))^2)) 
  
  result.perd.meancrps[i,1] <- mean(crps(f2(x), predy, predsig2))
  result.perd.meancrps[i,2] <- mean(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) 
  
  result.perd.comptime[i,1] <- toc.closed - tic.closed
  result.perd.comptime[i,2] <- toc.cokm - tic.cokm
}


py_run_file("Perdikaris.py")
result.perd.rmse <- cbind(result.perd.rmse, NARGP=unlist(py$l2error))
result.perd.meancrps <- cbind(result.perd.meancrps, NARGP=unlist(py$meancrps))
result.perd.comptime <- cbind(result.perd.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
### RMSE comparison ###
apply(result.perd.rmse, 2, mean) 
table(apply(result.perd.rmse, 1, which.min))
boxplot(result.perd.rmse)

### CRPS comparison ###
apply(result.perd.meancrps, 2, mean)
table(apply(result.perd.meancrps, 1, which.min))
boxplot(result.perd.meancrps)




