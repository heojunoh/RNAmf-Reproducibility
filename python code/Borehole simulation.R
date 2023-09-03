### Borehole Example ###
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
borehole <- function(xx)
{
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

boreholelow <- function(xx)
{ 
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 5 * Tu * (Hu-Hl) #+ (Tu * Kw) # Tu * Kw is added
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1.5+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

output.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  borehole(x[1:8])
} 

outputlow.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  boreholelow(x[1:8])
} 

### training data ###
n1 <- 80; n2 <- 40
d <- 8

rep <- 100
result.borehole.rmse <- matrix(NA, rep, 2)
# result.borehole.meanscore <- matrix(NA, rep, 4)
# result.borehole.medscore <- matrix(NA, rep, 4)
result.borehole.meancrps <- matrix(NA, rep, 2)
# result.borehole.medcrps <- matrix(NA, rep, 4)
result.borehole.comptime <- matrix(NA, rep, 2)
colnames(result.borehole.rmse) <- c("closed", "Cokriging")
# colnames(result.borehole.meanscore) <- c("closed", "Cokriging") # The larger, the better
# colnames(result.borehole.medscore) <- c("closed", "Cokriging") # The larger, the better
colnames(result.borehole.meancrps) <- c("closed", "Cokriging") # The smaller, the better
# colnames(result.borehole.medcrps) <- c("closed", "Cokriging") # The smaller, the better
colnames(result.borehole.comptime) <- c("closed", "Cokriging") # The smaller, the better

for(i in 1:rep) {
  set.seed(i)
  print(i)
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,outputlow.f)
  y2 <- apply(X2,1,output.f)

  
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
  
  ### compared to single fidelity ###
  # fit.GP2 <- GP(X2, y2)
  # pred2 <- pred.GP(fit.GP2, x)
  # 
  # ### direct fitting; not using closed form. f1(u) from (u, f1(u)) is random variable.
  # x1.mu <- rnorm(nrow(x), mean=pred.GP(fit.GP1, x)$mu, sd=sqrt(pred.GP(fit.GP1, x)$sig2)) 
  # xnew <- cbind(x, x1.mu) # Use mu of the input in the closed form
  # pred2new <- pred.GP(fit.GP2new, xnew) # not closed form
  
  # ### KOH method ###
  # y1d2 <- apply(X2,1,outputlow.f)
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
  # V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))+diag(eps,nrow(X1)+nrow(X2))
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
  
  
  ### RMSE ###
  # result.borehole.rmse[i,1] <- sqrt(mean((pred2$mu-apply(x,1,output.f))^2))/sqrt(sum((apply(x,1,output.f))^2)) # single fidelity
  result.borehole.rmse[i,1] <- sqrt(mean((predy-apply(x,1,output.f))^2)) # closed form
  # result.borehole.rmse[i,3] <- sqrt(mean((pred2new$mu-apply(x,1,output.f))^2))/sqrt(sum((apply(x,1,output.f))^2)) # not closed form
  result.borehole.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.f))^2)) # Cokm

  # result.borehole.meanscore[i,1] <- mean(score(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
  # result.borehole.meanscore[i,2] <- mean(score(apply(x,1,output.f), predy, predsig2)) # closed form
  # result.borehole.meanscore[i,3] <- mean(score(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
  # result.borehole.meanscore[i,4] <- mean(score(apply(x,1,output.f), pred.muficokm$mean, pred.muficokm$sig2)) # Cokm
  # 
  # result.borehole.medscore[i,1] <- median(score(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
  # result.borehole.medscore[i,2] <- median(score(apply(x,1,output.f), predy, predsig2)) # closed form
  # result.borehole.medscore[i,3] <- median(score(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
  # result.borehole.medscore[i,4] <- median(score(apply(x,1,output.f), pred.muficokm$mean, pred.muficokm$sig2)) # Cokm
  
  # result.borehole.meancrps[i,1] <- mean(crps(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
  result.borehole.meancrps[i,1] <- mean(crps(apply(x,1,output.f), predy, predsig2)) # closed form
  # result.borehole.meancrps[i,3] <- mean(crps(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
  result.borehole.meancrps[i,2] <- mean(crps(apply(x,1,output.f), pred.muficokm$mean, pred.muficokm$sig2)) # Cokm
  
  # result.borehole.medcrps[i,1] <- median(crps(apply(x,1,output.f), pred2$mu, pred2$sig2)) # single fidelity
  # result.borehole.medcrps[i,2] <- median(crps(apply(x,1,output.f), predy, predsig2)) # closed form
  # result.borehole.medcrps[i,3] <- median(crps(apply(x,1,output.f), pred2new$mu, pred2new$sig2)) # not closed form
  # result.borehole.medcrps[i,4] <- median(crps(apply(x,1,output.f), pred.muficokm$mean, pred.muficokm$sig2)) # Cokm
  
  result.borehole.comptime[i,1] <- toc.closed - tic.closed
  result.borehole.comptime[i,2] <- toc.cokm - tic.cokm
}

install.packages("reticulate")
library(reticulate)
py_run_file("Borehole.py")
result.borehole.rmse <- cbind(result.borehole.rmse, NARGP=unlist(py$l2error))
result.borehole.meancrps <- cbind(result.borehole.meancrps, NARGP=unlist(py$meancrps))
result.borehole.comptime <- cbind(result.borehole.comptime, NARGP=unlist(py$comptime))

par(mfrow=c(1,1))
#RMSE comparison#
apply(result.borehole.rmse, 2, mean) # 1.702758, 2.264872
table(apply(result.borehole.rmse, 1, which.min))
boxplot(result.borehole.rmse)

# #score comparison, The larger, the better
# apply(result.borehole.meanscore, 2, mean)
# table(apply(result.borehole.meanscore, 1, which.max))
# boxplot(result.borehole.meanscore)
# 
# #score comparison, The larger, the better
# apply(result.borehole.medscore, 2, mean)
# table(apply(result.borehole.medscore, 1, which.max))
# boxplot(result.borehole.medscore)

#CRPS comparison, The smaller, the better
apply(result.borehole.meancrps, 2, mean)
table(apply(result.borehole.meancrps, 1, which.min))
boxplot(result.borehole.meancrps)

# #CRPS comparison, The smaller, the better
# apply(result.borehole.medcrps, 2, mean)
# table(apply(result.borehole.medcrps, 1, which.min))
# boxplot(result.borehole.medcrps)



