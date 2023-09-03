### Currin Example ###
library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(RNAmf)

crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

### synthetic function ###
curretal88exp <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- 1 - exp(-1/(2*x2))
  fact2 <- 2300*x1^3 + 1900*x1^2 + 2092*x1 + 60
  fact3 <- 100*x1^3 + 500*x1^2 + 4*x1 + 20
  
  y <- fact1 * fact2/fact3
  return(y)
}

curretal88explc <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  maxarg <- max(c(0, x2-1/20))
  
  yh1 <- curretal88exp(c(x1+1/20, x2+1/20))
  yh2 <- curretal88exp(c(x1+1/20, maxarg))
  yh3 <- curretal88exp(c(x1-1/20, x2+1/20))
  yh4 <- curretal88exp(c(x1-1/20, maxarg))
  
  y <- (yh1 + yh2 + yh3 + yh4) / 4
  return(y)
}

### training data ###
n1 <- 20; n2 <- 10
d <- 2

rep <- 100
result.currin.rmse <- matrix(NA, rep, 2)
result.currin.meancrps <- matrix(NA, rep, 2)
result.currin.comptime <- matrix(NA, rep, 2)
colnames(result.currin.rmse) <- c("closed", "Cokriging")
colnames(result.currin.meancrps) <- c("closed", "Cokriging") # The smaller, the better
colnames(result.currin.comptime) <- c("closed", "Cokriging") # The smaller, the better


for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,curretal88explc)
  y2 <- apply(X2,1,curretal88exp)
  
  eps <- sqrt(.Machine$double.eps)
  
  ### test data ###
  x <- maximinLHS(1000, d)
  
  
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
  
  
  ### RMSE ###
  
  result.currin.rmse[i,1] <- sqrt(mean((predy-apply(x,1,curretal88exp))^2)) # closed form
  result.currin.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,curretal88exp))^2)) # Cokriging
  
  result.currin.meancrps[i,1] <- mean(crps(apply(x,1,curretal88exp), predy, predsig2)) # closed form
  result.currin.meancrps[i,2] <- mean(crps(apply(x,1,curretal88exp), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.currin.comptime[i,1] <- toc.closed - tic.closed
  result.currin.comptime[i,2] <- toc.cokm - tic.cokm
  
}

install.packages("reticulate")
library(reticulate)
py_run_file("Currin.py")
result.currin.rmse <- cbind(result.currin.rmse, NARGP=unlist(py$l2error))
result.currin.meancrps <- cbind(result.currin.meancrps, NARGP=unlist(py$meancrps))
result.currin.comptime <- cbind(result.currin.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.currin.rmse, 2, mean)
table(apply(result.currin.rmse, 1, which.min))
boxplot(result.currin.rmse)

#CRPS comparison, The smaller, the better
apply(result.currin.meancrps, 2, mean)
table(apply(result.currin.meancrps, 1, which.min))
boxplot(result.currin.meancrps)

