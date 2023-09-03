### nonlinear Example ###
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
result.nonlinear.rmse <- matrix(NA, rep, 2)
result.nonlinear.rmse2 <- matrix(NA, rep, 2)
# result.nonlinear.meanscore <- matrix(NA, rep, 2)
# result.nonlinear.medscore <- matrix(NA, rep, 2)
result.nonlinear.meancrps <- matrix(NA, rep, 2)
# result.nonlinear.medcrps <- matrix(NA, rep, 2)
result.nonlinear.comptime <- matrix(NA, rep, 2)
colnames(result.nonlinear.rmse) <- c("closed", "Cokriging")
# colnames(result.nonlinear.meanscore) <- c("closed", "Cokriging") # The larger, the better
# colnames(result.nonlinear.medscore) <- c("closed", "Cokriging") # The larger, the better
colnames(result.nonlinear.meancrps) <- c("closed", "Cokriging") # The smaller, the better
# colnames(result.nonlinear.medcrps) <- c("closed", "Cokriging") # The smaller, the better
colnames(result.nonlinear.comptime) <- c("closed", "Cokriging") # The smaller, the better

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
  
  
  # ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  # fit.GP1 <- GP(X1, y1, constant=TRUE)
  # 
  # ### model fitting using (x2, f1(x2)) ###
  # w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  # X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  # fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  
  ### test data ###
  x <- seq(0,1,length.out=1000)
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.closed <- predRNAmf(fit.closed, x)
  predy <- pred.closed$mu
  predsig2 <- pred.closed$sig2
  toc.closed <- proc.time()[3]
  
  # ### KOH method ###
  # fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  # pred.KOH2 <- pred.KOH(fit.KOH2, x)
  # mx1 <- pred.KOH2$mu
  # koh.var1 <- pred.KOH2$sig2


  ### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
  # w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
  # w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
  # xxnew <- cbind(x, w2.x)
  # pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
  # w1.x <- c(rep(NA, length(x)))
  # for(j in 1:length(x)){
  #   w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x[j])$mu, sd=sqrt(pred.GP(fit.GP1, x[j])$sig2)))
  # }

  # xxnew <- cbind(x, w1.x)
  # pred2new <- pred.GP(fit.GP2new, xxnew) # not closed form
  # 
  # ### prediction of original GP with single fidelity ###
  # fit.GP2 <- GP(X2, y2, constant=TRUE)
  # pred2 <- pred.GP(fit.GP2, x)

  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  ### RMSE ###
  # result.nonlinear.rmse[i,1] <- sqrt(sum((pred2$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # single fidelity
  # result.nonlinear.rmse[i,1] <- sqrt(sum((predy-f2(x))^2))/(sqrt(sum((f2(x))^2))) # closed form
  # result.nonlinear.rmse[i,3] <- sqrt(sum((pred2new$mu-f2(x))^2))/(sqrt(sum((f2(x))^2))) # not closed form
  # result.nonlinear.rmse[i,2] <- sqrt(sum((pred.muficokm$mean-f2(x))^2))/(sqrt(sum((f2(x))^2))) # Cokriging
  # result.nonlinear.rmse[i,5] <- sqrt(sum((mx1-f2(x))^2))/(sqrt(sum((f2(x))^2))) # KOH
  
  result.nonlinear.rmse[i,1] <- sqrt(mean((predy-f2(x))^2)) # closed form
  result.nonlinear.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-f2(x))^2)) # Cokriging
  
  # result.nonlinear.meanscore[i,1] <- mean(score(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  # result.nonlinear.meanscore[i,1] <- mean(score(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.meanscore[i,3] <- mean(score(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.meanscore[i,2] <- mean(score(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.meanscore[i,5] <- mean(score(f2(x), mx1, koh.var1)) # KOH

  # result.nonlinear.medscore[i,1] <- median(score(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  # result.nonlinear.medscore[i,1] <- median(score(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.medscore[i,3] <- median(score(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.medscore[i,2] <- median(score(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.medscore[i,5] <- median(score(f2(x), mx1, koh.var1)) # KOH
  
  # result.nonlinear.meancrps[i,1] <- mean(crps(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  result.nonlinear.meancrps[i,1] <- mean(crps(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.meancrps[i,3] <- mean(crps(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  result.nonlinear.meancrps[i,2] <- mean(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.meancrps[i,5] <- mean(crps(f2(x), mx1, koh.var1)) # KOH
  
  # result.nonlinear.medcrps[i,1] <- median(crps(f2(x), pred2$mu, pred2$sig2)) # single fidelity
  # result.nonlinear.medcrps[i,1] <- median(crps(f2(x), predy, predsig2)) # closed form
  # result.nonlinear.medcrps[i,3] <- median(crps(f2(x), pred2new$mu, pred2new$sig2)) # not closed form
  # result.nonlinear.medcrps[i,2] <- median(crps(f2(x), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.nonlinear.medcrps[i,5] <- median(crps(f2(x), mx1, koh.var1)) # KOH
  
  result.nonlinear.comptime[i,1] <- toc.closed - tic.closed
  result.nonlinear.comptime[i,2] <- toc.cokm - tic.cokm
}

install.packages("reticulate")
library(reticulate)
py_run_file("Perdikaris.py")
result.nonlinear.rmse <- cbind(result.nonlinear.rmse, NARGP=unlist(py$l2error))
result.nonlinear.meancrps <- cbind(result.nonlinear.meancrps, NARGP=unlist(py$meancrps))
result.nonlinear.comptime <- cbind(result.nonlinear.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.nonlinear.rmse, 2, mean) # 0.2825692, 0.4260575
table(apply(result.nonlinear.rmse, 1, which.min))
boxplot(result.nonlinear.rmse)

#score comparison, The larger, the better
# apply(result.nonlinear.meanscore, 2, mean)
# table(apply(result.nonlinear.meanscore, 1, which.max))
# boxplot(result.nonlinear.meanscore)
# 
# #score comparison, The larger, the better
# apply(result.nonlinear.medscore, 2, mean)
# table(apply(result.nonlinear.medscore, 1, which.max))
# boxplot(result.nonlinear.medscore)

#CRPS comparison, The smaller, the better
apply(result.nonlinear.meancrps, 2, mean)
table(apply(result.nonlinear.meancrps, 1, which.min))
boxplot(result.nonlinear.meancrps)

# #CRPS comparison, The smaller, the better
# apply(result.nonlinear.medcrps, 2, mean)
# table(apply(result.nonlinear.medcrps, 1, which.min))
# boxplot(result.nonlinear.medcrps)



