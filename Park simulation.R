### Park Example ###
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
park91a <- function(xx)
{
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
colnames(result.park.rmse) <- c("closed", "Cokm")
# result.park.meanscore <- matrix(NA, rep, 2)
# colnames(result.park.meanscore) <- c("closed", "Cokm")
result.park.meancrps <- matrix(NA, rep, 2)
colnames(result.park.meancrps) <- c("closed", "Cokm")
result.park.comptime <- matrix(NA, rep, 2)
colnames(result.park.comptime) <- c("closed", "Cokm")
meanf <- c(rep(0,100))


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
  
  
  ### model fitting for f1 ###
  eps <- sqrt(.Machine$double.eps)
  # fit.GP1 <- GP(X1, y1)
  # 
  # ### model fitting using (x2, f1(x2)) ###
  # w1.x2 <- pred.GP(fit.GP1, X2)$mu # can interpolate; nested
  # X2new <- cbind(X2, w1.x2) # combine (X2, f1(x2))
  # fit.GP2new <- GP(X2new, y2) # model fitting for f_M(X2, f1(x2))
  
  
  ### test data ###
  x <- maximinLHS(1000, d)
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.closed <- predRNAmf(fit.closed, x)
  predy <- pred.closed$mu
  predsig2 <- pred.closed$sig2
  toc.closed <- proc.time()[3]
  
  # ### compared to single fidelity ###
  # fit.GP2 <- GP(X2, y2)
  # pred2 <- pred.GP(fit.GP2, x)
  # 
  # ### direct fitting; not using closed form. f1(u) from (u, f1(u)) is random variable.
  # x1.mu <- rnorm(nrow(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
  # xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
  # pred2new <- pred.GP(fit.GP2new, xnew) # not closed form
  
  # ### KOH method ###
  # y1d2 <- apply(X2,1,park91blc)
  # 
  # ### estimating first order ###
  # fit.KOHGP1 <- KOHGP(X1, y1)
  # b1 <- 1/fit.KOHGP1$theta
  # sig2_1 <- fit.KOHGP1$tau2hat
  # 
  # ### estimating second order ###
  # # KOH(X2, y2, y1d2)
  # rho1 <- KOH(X2, y2, y1d2)$rho
  # b2 <- 1/KOH(X2, y2, y1d2)$theta
  # sig2_2 <- KOH(X2, y2, y1d2)$tau2hat
  # 
  # ### prediction of 2nd order KOH ###
  # tx1 <- cbind(rho1*sig2_1*covar.sep(x, X1, d=1/b1, g=eps), 
  #              rho1^2*sig2_1*covar.sep(x, X2, d=1/b1, g=eps) + sig2_2*covar.sep(x, X2, d=1/b2, g=eps))
  # 
  # V1 <- sig2_1*covar.sep(X1, d=1/b1, g=eps)
  # V12 <- rho1*sig2_1*covar.sep(X1, X2, d=1/b1, g=0)
  # V2 <- rho1^2*sig2_1*covar.sep(X2, d=1/b1, g=eps) + sig2_2*covar.sep(X2, d=1/b2, g=eps)
  # 
  # V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  # 
  # mx1 <- tx1 %*% solve(V_2) %*% c(y1, y2)
  # 
  # ### posterior variance ###
  # koh.var1 <- pmax(0, diag(sig2_2*covar.sep(x, d=1/b2, g=eps) + sig2_1*rho1^2*covar.sep(x, d=1/b1, g=eps) - tx1 %*% solve(V_2)%*%t(tx1)))  
  
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, #upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  meanf[i] <- sqrt(sum((apply(x,1,park91a))^2))
  
  ### RMSE ###
  # result.park.rmse[i,1] <- sqrt(mean((pred2$mu-apply(x,1,park91a))^2))/sqrt(sum((apply(x,1,park91a))^2)) # single fidelity
  result.park.rmse[i,1] <- sqrt(mean((predy-apply(x,1,park91a))^2)) # closed form
  # result.park.rmse[i,3] <- sqrt(mean((pred2new$mu-apply(x,1,park91a))^2))/sqrt(sum((apply(x,1,park91a))^2)) # not closed form
  result.park.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,park91a))^2)) # Cokm
  
  # result.park.meanscore[i,1] <- mean(score(apply(x,1,park91a), pred2$mu, pred2$sig2)) # single fidelity
  # result.park.meanscore[i,1] <- mean(score(apply(x,1,park91a), predy, predsig2)) # closed form
  # result.park.meanscore[i,3] <- mean(score(apply(x,1,park91a), pred2new$mu, pred2new$sig2)) # not closed form
  # result.park.meanscore[i,2] <- mean(score(apply(x,1,park91a), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  # result.park.meancrps[i,1] <- mean(crps(apply(x,1,park91a), pred2$mu, pred2$sig2)) # single fidelity
  result.park.meancrps[i,1] <- mean(crps(apply(x,1,park91a), predy, predsig2)) # closed form
  # result.park.meancrps[i,3] <- mean(crps(apply(x,1,park91a), pred2new$mu, pred2new$sig2)) # not closed form
  result.park.meancrps[i,2] <- mean(crps(apply(x,1,park91a), pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.park.comptime[i,1] <- toc.closed - tic.closed
  result.park.comptime[i,2] <- toc.cokm - tic.cokm
}

install.packages("reticulate")
library(reticulate)
py_run_file("Park.py")
result.park.rmse <- cbind(result.park.rmse, NARGP=unlist(py$l2error))
result.park.meancrps <- cbind(result.park.meancrps, NARGP=unlist(py$meancrps))
result.park.comptime <- cbind(result.park.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.park.rmse, 2, mean) # 0.05398435, 0.05866933
table(apply(result.park.rmse, 1, which.min))
boxplot(result.park.rmse)

#mean score#
apply(result.park.meanscore, 2, mean)
table(apply(result.park.meanscore, 1, which.min))
boxplot(result.park.meanscore)
c(sort(result.park.meanscore[,2])[6], apply(result.park.meanscore, 2, mean)[2], sort(result.park.meanscore[,2])[95])
c(sort(result.park.meanscore[,4])[6], apply(result.park.meanscore, 2, mean)[4], sort(result.park.meanscore[,4])[95])

#mean CRPS#
apply(result.park.meancrps, 2, mean)
table(apply(result.park.meancrps, 1, which.min))
boxplot(result.park.meancrps)
c(sort(result.park.meancrps[,2])[6], apply(result.park.meancrps, 2, mean)[2], sort(result.park.meancrps[,2])[95])
c(sort(result.park.meancrps[,4])[6], apply(result.park.meancrps, 2, mean)[4], sort(result.park.meancrps[,4])[95])



