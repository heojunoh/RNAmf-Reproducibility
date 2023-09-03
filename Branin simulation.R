### branin Example ###
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
branin <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  
  (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+ 10
}

braninm <- function(xx)
{ 
  x1 <- xx[1]
  x2 <- xx[2]
  
  10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 + (10-5/(4*pi))*cos((x1+2))+ 10) + 2*(x1-0.5) - 3*(3*x2-1) - 1
}

braninl <- function(xx)
{ x1 <- xx[1]
  x2 <- xx[2]
  
  10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 + (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) + 2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
}

output.branin <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  branin(x[1:2])
} 

output.braninl <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  braninl(x[1:2])
} 

output.braninm <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  braninm(x[1:2])
} 


### training data ###
n1 <- 20; n2 <- 15; n3 <- 10

rep <- 100
result.branin.rmse <- matrix(NA, rep, 2)
# result.branin.meanscore <- matrix(NA, rep, 2)
# result.branin.medscore <- matrix(NA, rep, 2)
result.branin.meancrps <- matrix(NA, rep, 2)
# result.branin.medcrps <- matrix(NA, rep, 2)
result.branin.comptime <- matrix(NA, rep, 2)
colnames(result.branin.rmse) <- c("closed", "Cokriging")
# colnames(result.branin.meanscore) <- c("closed", "Cokriging") # The larger, the better
# colnames(result.branin.medscore) <- c("closed", "Cokriging") # The larger, the better
colnames(result.branin.meancrps) <- c("closed", "Cokriging") # The smaller, the better
# colnames(result.branin.medcrps) <- c("closed", "Cokriging") # The smaller, the better
colnames(result.branin.comptime) <- c("closed", "Cokriging") # The smaller, the better
meanf <- c(rep(0,100))

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
  
  y1 <- apply(X1,1,output.braninl)
  y2 <- apply(X2,1,output.braninm)
  y3 <- apply(X3,1,output.branin)
  
  
  # ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  # fit.GP1 <- GP(X1, y1, constant=TRUE)
  # 
  # ### model fitting using (x2, f1(x2)) ###
  # w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  # X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  # fit.GP2new <- GP(X2new, y2, constant=TRUE) # model fitting for f_M(X2, f1(x2))
  # 
  # ### model fitting using (x3, f2(x3, f1(x3))) ###
  # w1.x3 <- pred.GP(fit.GP1, X3)$mu # can interpolate; nested
  # w2.x3 <- pred.GP(fit.GP2new, cbind(X3, w1.x3))$mu # can interpolate; nested
  # X3new <- cbind(X3, w2.x3) # combine (X3, f2(x3, f1(x3)))
  # fit.GP3new <- GP(X3new, y3, constant=TRUE) # model fitting for f_H(X3, f2(x3, f1(x3)))
  
  
  ### test data ###
  x <- maximinLHS(1000, 2)
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  pred.closed <- predRNAmf2(fit.closed, x)
  predy <- pred.closed$mu
  predsig2 <- pred.closed$sig2
  toc.closed <- proc.time()[3]
  
  ### KOH method ###
  # fit.KOH3 <- fit.KOH(X1, X2, X3, y1, y2, y3)
  # pred.KOH3 <- pred.KOH(fit.KOH3, x)
  # mx2 <- pred.KOH3$mu
  # koh.var2 <- pred.KOH3$sig2


  ### direct fitting; not using closed form. f1(u) and f_M(u) from (u, f_M(u, f1(u))) are random variables.
  # w1.x <- rnorm(length(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) # sample f1(x)
  # w2.x <- rnorm(length(x), mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu, sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2))
  # xxnew <- cbind(x, w2.x)
  # pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form
  # w1.x <- c(rep(NA, nrow(x)))
  # w2.x <- c(rep(NA, nrow(x)))
  # for(j in 1:nrow(x)){
  #   w1.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP1, x)$mu[j], sd=sqrt(pred.GP(fit.GP1, x)$sig2[j])))
  # }
  # for(j in 1:nrow(x)){
  #   w2.x[j] <- mean(rnorm(10000, mean=pred.GP(fit.GP2new, cbind(x, w1.x))$mu[j], sd=sqrt(pred.GP(fit.GP2new, cbind(x, w1.x))$sig2[j])))
  # }
  # 
  # xxnew <- cbind(x, w2.x)
  # pred3new <- pred.GP(fit.GP3new, xxnew) # not closed form

  # ### prediction of original GP with single fidelity ###
  # fit.GP3 <- GP(X3, y3, constant=TRUE)
  # pred3 <- pred.GP(fit.GP3, x)

  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.5,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]

  meanf[i] <- (sqrt(sum((apply(x,1,output.branin))^2)))
  
  ### RMSE ###
  # result.branin.rmse[i,1] <- sqrt(sum((pred3$mu-apply(x,1,output.branin))^2))/(sqrt(sum((apply(x,1,output.branin))^2))) # single fidelity
  result.branin.rmse[i,1] <- sqrt(mean((predy-apply(x,1,output.branin))^2)) # closed form
  # result.branin.rmse[i,3] <- sqrt(sum((pred3new$mu-apply(x,1,output.branin))^2))/(sqrt(sum((apply(x,1,output.branin))^2))) # not closed form
  result.branin.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.branin))^2)) # Cokriging
  # result.branin.rmse[i,5] <- sqrt(sum((mx2-apply(x,1,output.branin))^2))/(sqrt(sum((apply(x,1,output.branin))^2))) # KOH

  # result.branin.meanscore[i,1] <- mean(score(apply(x,1,output.branin), pred3$mu, pred3$sig2)) # single fidelity
  # result.branin.meanscore[i,1] <- mean(score(apply(x,1,output.branin), predy, predsig2)) # closed form
  # result.branin.meanscore[i,3] <- mean(score(apply(x,1,output.branin), pred3new$mu, pred3new$sig2)) # not closed form
  # result.branin.meanscore[i,2] <- mean(score(apply(x,1,output.branin), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.branin.meanscore[i,5] <- mean(score(apply(x,1,output.branin), mx2, koh.var2)) # KOH
  
  # result.branin.medscore[i,1] <- median(score(apply(x,1,output.branin), pred3$mu, pred3$sig2)) # single fidelity
  # result.branin.medscore[i,1] <- median(score(apply(x,1,output.branin), predy, predsig2)) # closed form
  # result.branin.medscore[i,3] <- median(score(apply(x,1,output.branin), pred3new$mu, pred3new$sig2)) # not closed form
  # result.branin.medscore[i,2] <- median(score(apply(x,1,output.branin), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.branin.medscore[i,5] <- median(score(apply(x,1,output.branin), mx2, koh.var2)) # KOH
  
  # result.branin.meancrps[i,1] <- mean(crps(apply(x,1,output.branin), pred3$mu, pred3$sig2)) # single fidelity
  result.branin.meancrps[i,1] <- mean(crps(apply(x,1,output.branin), predy, predsig2)) # closed form
  # result.branin.meancrps[i,3] <- mean(crps(apply(x,1,output.branin), pred3new$mu, pred3new$sig2)) # not closed form
  result.branin.meancrps[i,2] <- mean(crps(apply(x,1,output.branin), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.branin.meancrps[i,5] <- mean(crps(apply(x,1,output.branin), mx2, koh.var2)) # KOH
  
  # result.branin.medcrps[i,1] <- median(crps(apply(x,1,output.branin), pred3$mu, pred3$sig2)) # single fidelity
  # result.branin.medcrps[i,1] <- median(crps(apply(x,1,output.branin), predy, predsig2)) # closed form
  # result.branin.medcrps[i,3] <- median(crps(apply(x,1,output.branin), pred3new$mu, pred3new$sig2)) # not closed form
  # result.branin.medcrps[i,2] <- median(crps(apply(x,1,output.branin), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.branin.medcrps[i,5] <- median(crps(apply(x,1,output.branin), mx2, koh.var2)) # KOH
  
  result.branin.comptime[i,1] <- toc.closed - tic.closed
  result.branin.comptime[i,2] <- toc.cokm - tic.cokm
}

install.packages("reticulate")
library(reticulate)
py_run_file("Branin.py")
result.branin.rmse <- cbind(result.branin.rmse, NARGP=unlist(py$l2error))
result.branin.meancrps <- cbind(result.branin.meancrps, NARGP=unlist(py$meancrps))
result.branin.comptime <- cbind(result.branin.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.branin.rmse, 2, mean) # 30.93884, 30.97887
table(apply(result.branin.rmse, 1, which.min))
boxplot(result.branin.rmse)

# #score comparison, The larger, the better
# apply(result.branin.meanscore, 2, mean)
# table(apply(result.branin.meanscore, 1, which.max))
# boxplot(result.branin.meanscore)
# 
# #score comparison, The larger, the better
# apply(result.branin.medscore, 2, mean)
# table(apply(result.branin.medscore, 1, which.max))
# boxplot(result.branin.medscore)

#CRPS comparison, The smaller, the better
apply(result.branin.meancrps, 2, mean)
table(apply(result.branin.meancrps, 1, which.min))
boxplot(result.branin.meancrps)

# #CRPS comparison, The smaller, the better
# apply(result.branin.medcrps, 2, mean)
# table(apply(result.branin.medcrps, 1, which.min))
# boxplot(result.branin.medcrps)



