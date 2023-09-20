costmatc3 <- list(NA)
rmsematc3 <- list(NA)
crpsmatc3 <- list(NA)
time.each3 <- rep(0,10)
cost <- 3
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

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  print(kk)
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  ### test data ###
  x <- maximinLHS(1000, d)
  
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.RNAmf, x)$mu
  predsig2 <- predRNAmf(fit.RNAmf, x)$sig2
  
  ### RMSE ###
  print(sqrt(mean((predy-apply(x,1,park91a))^2))) # closed form
  
  park.cost <- 0
  park.error <- sqrt(mean((predy-apply(x,1,park91a))^2))
  park.crps <- mean(crps(apply(x,1,park91a), predy, predsig2))
  
  Iselect <- ALMC_two_level(x, fit.RNAmf, 100, c(1,cost), list(park91alc, park91a), parallel=TRUE, ncore=10)
  
  
  #################
  ### Add point ###
  #################
  while(park.cost[length(park.cost)] < 30){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, x)$mu
    predsig2 <- predRNAmf(Iselect$fit, x)$sig2
    
    ### RMSE ###
    park.error <- c(park.error, sqrt(mean((predy-apply(x,1,park91a))^2))) # RMSE
    park.crps <- c(park.crps, mean(crps(apply(x,1,park91a), predy, predsig2))) # CRPS
    if(Iselect$chosen$level == 1){
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+1
    }else{
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+(1+cost)
    }
    print(park.cost[length(park.cost)])
    print(park.error[length(park.error)])
    
    if(park.cost[length(park.cost)] >= 30){break}
    
    ### update the next point ###
    Iselect <- ALMC_two_level(x, Iselect$fit, 100, c(1,cost), list(park91alc, park91a), parallel=TRUE, ncore=10)
    # save.image("Park AL3 1,3.RData")
  }
  
  ### Save results ###
  costmatc3[[kk]] <- park.cost
  rmsematc3[[kk]] <- park.error
  crpsmatc3[[kk]] <- park.crps
  
  time.each3[kk] <- proc.time()[3]- time.start
  # save.image("Park AL3 1,3.RData")
}
costmatc3
rmsematc3
crpsmatc3
time.each3
