costmatc3 <- list(NA)
rmsematc3 <- list(NA)
crpsmatc3 <- list(NA)
time.each3 <- rep(0,10)
cost <- c(1,3)
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

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  cat("Perdikaris - ALMC: ", kk, "/ 10\n")
  
  X <- NestedX(c(n1, n2), 1)
  X1 <- X[[1]]
  X2 <- X[[2]]
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  ### test data ###
  x <- seq(0,1,length.out=n.test)
  
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predict(fit.RNAmf, x)$mu
  predsig2 <- predict(fit.RNAmf, x)$sig2
  
  ### RMSE ###
  # sqrt(mean((predy-f2(x))^2)) # RNAmf
  
  perd.cost <- 0
  perd.error <- sqrt(mean((predy-f2(x))^2))
  perd.crps <- mean(crps(f2(x), predy, predsig2))
  
  almc.RNAmf <- ALMC_RNAmf(Xref = x, Xcand = x, fit.RNAmf, cost = cost,
                           optim = TRUE, parallel = parallel, ncore = ncore)
  
  if (almc.RNAmf$chosen$level == 1) {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), almc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, f1(almc.RNAmf$chosen$Xnext)),
                                      matrix(t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,1]),
                                      fit.RNAmf$fit2$y,
                                      kernel="sqex", constant=TRUE)
  } else {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), almc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, f1(almc.RNAmf$chosen$Xnext)),
                                      rbind(matrix(t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,1]), almc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit2$y, f2(almc.RNAmf$chosen$Xnext)),
                                      kernel="sqex", constant=TRUE)
  }
  
  #################
  ### Add point ###
  #################
  while(perd.cost[length(perd.cost)] < 50){ # if total cost is less than the budget
    ### RNAmf ###
    predy <- predict(fit.RNAmf.next, x)$mu
    predsig2 <- predict(fit.RNAmf.next, x)$sig2
    
    ### RMSE ###
    perd.error <- c(perd.error, sqrt(mean((predy-f2(x))^2))) # RMSE
    perd.crps <- c(perd.crps, mean(crps(f2(x), predy, predsig2))) # CRPS
    if(almc.RNAmf$chosen$level == 1){
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+cost[1]
    }else{
      perd.cost[length(perd.cost)+1] <- perd.cost[length(perd.cost)]+sum(cost)
    }
    # print(perd.cost[length(perd.cost)])
    # print(perd.error[length(perd.error)])
    
    # if(perd.cost[length(perd.cost)] >= 50){break}

    ### update the next point ###
    almc.RNAmf <- ALMC_RNAmf(Xref = x, Xcand = x, fit.RNAmf.next, cost = cost,
                             optim = TRUE, parallel = parallel, ncore = ncore)
    
    if (almc.RNAmf$chosen$level == 1) {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), almc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, f1(almc.RNAmf$chosen$Xnext)),
                                        matrix(t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,1]),
                                        fit.RNAmf.next$fit2$y,
                                        kernel="sqex", constant=TRUE)
    } else {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), almc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, f1(almc.RNAmf$chosen$Xnext)),
                                        rbind(matrix(t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,1]), almc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit2$y, f2(almc.RNAmf$chosen$Xnext)),
                                        kernel="sqex", constant=TRUE)
    }
    
    # save.image("/Perd AL3 1,3.RData")
  }
  
  ### Save results ###
  costmatc3[[kk]] <- perd.cost
  rmsematc3[[kk]] <- perd.error
  crpsmatc3[[kk]] <- perd.crps
  
  time.each3[kk] <- proc.time()[3]- time.start
  # save.image("Perd ALMC.RData")
}
costmatc3
rmsematc3
crpsmatc3
time.each3
