costmatk <- list(NA)
rmsematk <- list(NA)
crpsmatk <- list(NA)
cost <- 3
### Perdikaris function ###
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
d <- 1
### test data and number of candidate points ###
n.test <- 1000

KOHGP <- function(X, y, g=eps, center=TRUE){ # For 1st level's hyperparameter
  
  n <- length(y)
  # if(center) y <- scale(y, scale=FALSE)
  
  nlsep <- function(par, X, Y) 
  {
    n <- length(Y)
    theta <- par # lengthscale
    K <- covar.sep(X, d=theta, g=g)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(drop(-ll))
  }
  
  outg <- optim(rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X)), nlsep, method="L-BFGS-B", 
                lower=c(rep(0.001*sqrt(ncol(X)), ncol(X))), upper=c(rep(1000*sqrt(ncol(X)), ncol(X))), X=X, Y=y) 
  
  K <- covar.sep(X, d=outg$par, g=g)
  Ki <- solve(K+diag(g,n))
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

KOH <- function(X, Y2, Y1, g=eps, center=TRUE){ # For after 1st level's hyperparameter
  
  n <- length(Y2)
  # if(center) y <- scale(y, scale=FALSE)
  
  nlsep2 <- function(par, X, Y2, Y1) 
  {
    n <- length(Y2)
    theta <- par # lengthscale and rho
    Y <- Y2-theta[ncol(X)+1]*Y1 # y2-rho*y1
    K <- covar.sep(X, d=theta[1:ncol(X)], g=g)
    Ki <- solve(K+diag(g,n))
    ldetK <- determinant(K, logarithm=TRUE)$modulus
    ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
    return(drop(-ll))
  }
  
  outg <- optim(c(rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X)), 0.1), # For lengthscale and rho #
                nlsep2, method="L-BFGS-B", 
                lower=c(rep(0.001*sqrt(ncol(X)), ncol(X)+1)), 
                upper=c(rep(1000*sqrt(ncol(X)), ncol(X)+1)), 
                X=X, Y2=Y2, Y1=Y1) 
  
  K <- covar.sep(X, d=outg$par[1:ncol(X)], g=g)
  Ki <- solve(K+diag(g,n))
  rho <- outg$par[ncol(X)+1]
  y <- Y2-rho*Y1
  tau2hat <- drop(t(y) %*% Ki %*% y / n)
  
  return(list(theta = outg$par[1:ncol(X)], rho = rho, g=g, Ki=Ki, X = X, y = y, tau2hat=tau2hat))
}

IMSPEKOH1 <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  # level; level of fidelity
  X1 <- fit$X1
  X2 <- fit$X2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  if(level==1){ ### level 1 ###
    X1 <- rbind(fit$X1, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:d], g=g),
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[(d+1):(2*d)], g=g))
    
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:d], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:d], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[(d+1):(2*d)], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[(d+1):(2*d)], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:d], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
    
  }else if(level==2){ ### level 2 ###
    X2 <- rbind(fit$X2, newx)
    
    ### update sig2
    tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:d], g=g),
                 rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[(d+1):(2*d)], g=g))
    
    V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:d], g=g)
    V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:d], g=0)
    V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[(d+1):(2*d)], g=g)
    
    V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
    
    koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[(d+1):(2*d)], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:d], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  }else{
    stop("level should be 1 or 2.")
  }
  
  return(IMSPE=mean(koh.var1))
}

IMSPEKOHselect1 <- function(x, newx, fit, level){
  ### This is updating KOH when the one data point is added ###
  # x; usually grid
  # newx; new point which will be added, should be n=1
  # fit; fitted KOH model
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  # b <- fit$b
  # rho <- fit$rho
  # tau2hat <- fit$tau2hat
  # g <- fit$g
  
  if(level==1){ ### level 1 ###
    ### Generate output
    y.select <- f1(newx)
    
    X1 <- rbind(fit$X1, newx)
    Y1 <- c(fit$Y1, y.select)
    
  }else if(level==2){ ### level 2 ###
    ### Generate output
    y.select <- f2(newx)
    
    X2 <- rbind(fit$X2, newx)
    Y2 <- c(fit$Y2, y.select)
    
  }else{
    stop("level should be 1 or 2.")
  }
  
  fit <- fit.KOH(X1, X2, Y1, Y2)
  
  return(fit=fit)
}

fit.KOH <- function(X1, X2, Y1, Y2, g=eps){ # need to change function for another example
  
  ### KOH method ###
  Y1d2 <- f1(X2)
  
  ### estimating first order ###
  fit.KOHGP1 <- KOHGP(X1, Y1)
  b1 <- 1/fit.KOHGP1$theta
  sig2_1 <- fit.KOHGP1$tau2hat
  
  ### estimating second order ###
  # KOH(X2, Y2, Y1d2)
  rho1 <- KOH(X2, Y2, Y1d2)$rho
  b2 <- 1/KOH(X2, Y2, Y1d2)$theta
  sig2_2 <- KOH(X2, Y2, Y1d2)$tau2hat
  
  return(list(b=c(b1, b2), rho=rho1, tau2hat=c(sig2_1, sig2_2), g=g, X1=X1, X2=X2, Y1=Y1, Y2=Y2))
}

pred.KOH <- function(fit, x){ # need to change function for another example
  
  X1 <- fit$X1
  X2 <- fit$X2
  Y1 <- fit$Y1
  Y2 <- fit$Y2
  
  b <- fit$b
  rho <- fit$rho
  tau2hat <- fit$tau2hat
  g <- fit$g
  
  
  ### prediction of 2nd order KOH ###
  tx1 <- cbind(rho*tau2hat[1]*covar.sep(x, X1, d=1/b[1:d], g=g), 
               rho^2*tau2hat[1]*covar.sep(x, X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(x, X2, d=1/b[(d+1):(2*d)], g=g))
  
  V1 <- tau2hat[1]*covar.sep(X1, d=1/b[1:d], g=g)
  V12 <- rho*tau2hat[1]*covar.sep(X1, X2, d=1/b[1:d], g=0)
  V2 <- rho^2*tau2hat[1]*covar.sep(X2, d=1/b[1:d], g=g) + tau2hat[2]*covar.sep(X2, d=1/b[(d+1):(2*d)], g=g)
  
  V_2 <- rbind(cbind(V1, V12), cbind(t(V12), V2))
  
  mx1 <- tx1 %*% solve(V_2) %*% c(Y1, Y2)
  
  ### posterior variance ###
  koh.var1 <- pmax(0, diag(tau2hat[2]*covar.sep(as.matrix(x), d=1/b[(d+1):(2*d)], g=g) + tau2hat[1]*rho^2*covar.sep(as.matrix(x), d=1/b[1:d], g=g) - tx1 %*% solve(V_2)%*%t(tx1)))
  
  return(list(mu=mx1, sig2=koh.var1))
}


for(kk in 1:10){
  set.seed(kk)
  cat("Perdikaris - MRSUR: ", kk, "/ 10\n")
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  ### test data ###
  x <- seq(0,1,length.out=n.test)
  
  ### KOH method ###
  fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  pred.KOH2 <- pred.KOH(fit.KOH2, x)
  mx1 <- pred.KOH2$mu
  koh.var1 <- pred.KOH2$sig2
  
  ### RMSE ###
  # sqrt(mean((mx1-f2(x))^2)) # KOH
  
  
  ### IMSPE ###
  g <- seq(0,1,length.out=n.test) # more than 1-dim, use expand.grid()
  Icurrent <- mean(koh.var1) # current IMSPE
  Icurrent
  
  ### Add 1 points and calculate IMSPE ###
  IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
  IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
  
  for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
    IcandKOH1[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=1)
  }
  for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
    IcandKOH2[i] <- IMSPEKOH1(g, g[i], fit.KOH2, level=2)
  }
  
  mrsur <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])

  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(mrsur/c(1,cost))
  chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2)[,chosen[1,1]])
  
  perd.cost <- 0
  perd.error <- sqrt(mean((mx1-f2(x))^2))
  perd.crps <- mean(crps(f2(x), mx1, koh.var1))
  
  Iselect <- IMSPEKOHselect1(g, g[chosen[nrow(chosen),2]], fit.KOH2, level=chosen[nrow(chosen),1]) 
  
  
  #################
  ### Add point ###
  #################
  while(perd.cost[length(perd.cost)] < 50){ # if total cost is less than the budget
    
    ### predictive ###
    pred.KOH2 <- pred.KOH(Iselect, x)
    mx1 <- pred.KOH2$mu
    koh.var1 <- pred.KOH2$sig2
    
    ### RMSE ###  
    perd.error <- c(perd.error, sqrt(mean((mx1-f2(x))^2))) # RMSE
    perd.crps <- c(perd.crps, mean(crps(f2(x), mx1, koh.var1))) # CRPS
    if(chosen[nrow(chosen),1] == 1){
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+1
    }else{
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+cost
    }
    
    #############
    ### IMSPE ###
    #############
    g <- seq(0,1,length.out=n.test) # more than 1-dim, use expand.grid()
    Icurrent <- mean(koh.var1)
    
    ### Add 1 points and calculate IMSPE ###
    IcandKOH1 <- c(rep(0, length(g))) # IMSPE candidates
    IcandKOH2 <- c(rep(0, length(g))) # IMSPE candidates
    
    for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
        IcandKOH1[i] <- IMSPEKOH1(g, g[i], Iselect, level=1)
      }
    }
    for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
        IcandKOH2[i] <- IMSPEKOH1(g, g[i], Iselect, level=2)
      }
    }
    
    if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
    if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
    
    mrsur <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
   
    chosen <- rbind(chosen, c(which.max(mrsur/c(1,cost)), which.min(cbind(IcandKOH1, IcandKOH2)[,which.max(mrsur/c(1,cost))])))
    Iselect <- IMSPEKOHselect1(g, g[chosen[nrow(chosen),2]], Iselect, level=chosen[nrow(chosen),1])
    
    if(perd.cost[length(perd.cost)] >= 50){break}
    
  }
  
  ### Save results ###
  costmatk[[kk]] <- perd.cost
  rmsematk[[kk]] <- perd.error
  crpsmatk[[kk]] <- perd.crps
}
costmatk
rmsematk
crpsmatk
