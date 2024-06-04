costmatk <- list(NA)
rmsematk <- list(NA)
crpsmatk <- list(NA)
cost <- 3

### Park function ###
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
    y.select <- park91alc(newx)
    
    X1 <- rbind(fit$X1, newx)
    Y1 <- c(fit$Y1, y.select)
    
  }else if(level==2){ ### level 2 ###
    ### Generate output
    y.select <- park91a(newx)
    
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
  Y1d2 <- apply(X2,1,park91alc) 
  
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
  cat("Park - MRSUR: ", kk, "/ 10\n")
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  ### test data ###
  x <- maximinLHS(n.test, d)
  
  ### KOH method ###
  fit.KOH2 <- fit.KOH(X1, X2, y1, y2)
  pred.KOH2 <- pred.KOH(fit.KOH2, x)
  mx1 <- pred.KOH2$mu
  koh.var1 <- pred.KOH2$sig2
  
  ### RMSE ###
  # sqrt(mean((mx1-apply(x,1,park91a))^2)) # KOH
  
  
  ### IMSPE ###
  Icurrent <- mean(koh.var1) # current IMSPE
  
  ### Add 1 points and calculate IMSPE ###
  IcandKOH1 <- c(rep(0, nrow(x))) # IMSPE candidates
  IcandKOH2 <- c(rep(0, nrow(x))) # IMSPE candidates
  
  for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
    IcandKOH1[i] <- IMSPEKOH1(x, x[i,], fit.KOH2, level=1)
  }
  for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
    IcandKOH2[i] <- IMSPEKOH1(x, x[i,], fit.KOH2, level=2)
  }

  mrsur <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
 
  chosen <- matrix(0, ncol=2)
  chosen[1,1] <- which.max(mrsur/c(1,cost))
  chosen[1,2] <- which.min(cbind(IcandKOH1, IcandKOH2)[,chosen[1,1]])
  

  park.cost <- 0
  park.error <- sqrt(mean((mx1-apply(x,1,park91a))^2))
  park.crps <- mean(crps(apply(x,1,park91a), mx1, koh.var1))
  
  Iselect <- IMSPEKOHselect1(x, x[chosen[nrow(chosen),2],], fit.KOH2, level=chosen[nrow(chosen),1]) 
  
  
  #################
  ### Add point ###
  #################
  while(park.cost[length(park.cost)] < 30){ # if total cost is less than the budget
    
    ### predictive ###
    pred.KOH2 <- pred.KOH(Iselect, x)
    mx1 <- pred.KOH2$mu
    koh.var1 <- pred.KOH2$sig2
    
    ### RMSE ###  
    park.error <- c(park.error, sqrt(mean((mx1-apply(x,1,park91a))^2))) # RMSE
    park.crps <- c(park.crps, mean(crps(apply(x,1,park91a), mx1, koh.var1))) # CRPS
    if(chosen[nrow(chosen),1] == 1){
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+1
    }else{
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+cost
    }
    
    #############
    ### IMSPE ###
    #############
    Icurrent <- mean(koh.var1)
    
    ### Add 1 points and calculate IMSPE ###
    IcandKOH1 <- c(rep(0, nrow(x))) # IMSPE candidates
    IcandKOH2 <- c(rep(0, nrow(x))) # IMSPE candidates
    
    for(i in 1:length(IcandKOH1)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH1[i] <- 0}else{
        IcandKOH1[i] <- IMSPEKOH1(x, x[i,], Iselect, level=1)
      }
    }
    for(i in 1:length(IcandKOH2)){ # no true, no need to fit just pred
      if(any(chosen[,2]==i)){IcandKOH2[i] <- 0}else{
        IcandKOH2[i] <- IMSPEKOH1(x, x[i,], Iselect, level=2)
      }
    }
    
    if(any(IcandKOH1==0)){IcandKOH1[which(IcandKOH1==0)] <-  max(IcandKOH1)}
    if(any(IcandKOH2==0)){IcandKOH2[which(IcandKOH2==0)] <-  max(IcandKOH2)}
    
    mrsur <- c(Icurrent - IcandKOH1[which.min(IcandKOH1)], Icurrent - IcandKOH2[which.min(IcandKOH2)])
    
    chosen <- rbind(chosen, c(which.max(mrsur/c(1,cost)), which.min(cbind(IcandKOH1, IcandKOH2)[,which.max(mrsur/c(1,cost))])))
    Iselect <- IMSPEKOHselect1(x, x[chosen[nrow(chosen),2],], Iselect, level=chosen[nrow(chosen),1])
    
    if(park.cost[length(park.cost)] >=30){break}
    
  }

  ### Save results ###
  costmatk[[kk]] <- park.cost
  rmsematk[[kk]] <- park.error
  crpsmatk[[kk]] <- park.crps
}
costmatk
rmsematk
crpsmatk
