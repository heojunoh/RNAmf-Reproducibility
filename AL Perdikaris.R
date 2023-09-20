eps <- sqrt(.Machine$double.eps)

costmatc <- list(NA)
rmsematc <- list(NA)
crpsmatc <- list(NA)
time.each <- rep(0,10)
cost <- 3
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

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  print(kk)
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)

  NestDesign <- NestedDesignBuild(design = list(X1,X2))

  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)

  y1 <- f1(X1)
  y2 <- f2(X2)

  ### test data ###
  x <- seq(0,1,length.out=1000)


  ### RNAmf ###
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.RNAmf, x)$mu
  predsig2 <- predRNAmf(fit.RNAmf, x)$sig2

  ### RMSE ###
  sqrt(mean((predy-f2(x))^2)) # RNAmf
  
  perd.cost <- 0
  perd.error <- sqrt(mean((predy-f2(x))^2))
  perd.crps <- mean(crps(f2(x), predy, predsig2))
  
  Iselect <- ALM_two_level(fit.RNAmf, c(1,cost), list(f1, f2), parallel=TRUE, ncore=10)


  #################
  ### Add point ###
  #################
  while(perd.cost[length(perd.cost)] < 43){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2

    ### RMSE ###
    perd.error <- c(perd.error, sqrt(mean((predy-f2(x))^2))) # RMSE
    perd.crps <- c(perd.crps, mean(crps(f2(x), predy, predsig2))) # CRPS
    if(Iselect$chosen$level == 1){
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+1
    }else{
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+(1+cost)
    }
    print(perd.cost[length(perd.cost)])
    print(perd.error[length(perd.error)])
    
    if(perd.cost[length(perd.cost)] >= 43){break}
    
    ### update the next point ###
    Iselect <- ALM_two_level(Iselect$fit, c(1,cost), list(f1, f2), parallel=TRUE, ncore=10)
    # save.image("/Perd AL 1,3.RData")
  }

  ### Save results ###
  costmatc[[kk]] <- perd.cost
  rmsematc[[kk]] <- perd.error
  crpsmatc[[kk]] <- perd.crps
  
  time.each[kk] <- proc.time()[3]- time.start
  # save.image("/Perd AL 1,3.RData")
}
costmatc
rmsematc
crpsmatc
time.each
