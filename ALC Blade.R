costmatc2 <- list(NA)
rmsematc2 <- list(NA)
crpsmatc2 <- list(NA)
time.each2 <- rep(0,10)
cost <- c(2.25, 6.85)

# Blade simulation
blade1 <- function(xx){
  d1 <- data.frame(xx*0.5+0.25, rep(0.05, nrow(xx))) # scale X to [-1,1]
  write.csv(xx, "Rmatlab_files/generate_text/temp_to_X.txt", row.names=F) # change path
  write.csv(d1, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F) # change path
  run_matlab_script("Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,  # change path
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",") # change path
  y <- d2$V4
  y
}

blade2 <- function(xx){
  d1 <- data.frame(xx*0.5+0.25, rep(0.0125, nrow(xx))) # scale X to [-1,1]
  write.csv(xx, "Rmatlab_files/generate_text/temp_to_X.txt", row.names=F) # change path
  write.csv(d1, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F) # change path
  run_matlab_script("Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,  # change path
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",") # change path
  y <- d2$V4
  y
}

### training data ###
n1 <- 20; n2 <- 10
d <- 2
### test data and number of candidate points ###
n.test <- 100

### Test data ###
set.seed(1)
X.test <- maximinLHS(n.test, d)
y.test <- blade2(X.test)

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  cat("Blade - ALC: ", kk, "/ 10\n")
  
  ### Generate Input ###
  X <- NestedX(c(n1, n2), d)
  X1 <- X[[1]]
  X2 <- X[[2]]
  
  ### Y1 ###
  y1 <- blade1(X1)
  
  ### Y2 ###
  y2 <- blade2(X2)
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predict(fit.RNAmf, X.test)$mu
  predsig2 <- predict(fit.RNAmf, X.test)$sig2
  
  ### RMSE ###
  # print(sqrt(mean((predy-y.test)^2))) 
  
  blade.cost <- 0
  blade.error <- sqrt(mean((predy-y.test)^2))
  blade.crps <- mean(crps(y.test, predy, predsig2))
    
  alc.RNAmf <- ALC_RNAmf(Xref = X.test, Xcand = X.test, fit.RNAmf, cost = cost, 
                         optim = TRUE, parallel = parallel, ncore = ncore)
  if (alc.RNAmf$chosen$level == 1) {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), alc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, blade1(alc.RNAmf$chosen$Xnext)),
                                      t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)],
                                      fit.RNAmf$fit2$y,
                                      kernel="sqex", constant=TRUE)
  } else {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), alc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, blade1(alc.RNAmf$chosen$Xnext)),
                                      rbind(t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)], alc.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit2$y, blade2(alc.RNAmf$chosen$Xnext)),
                                      kernel="sqex", constant=TRUE)
  }
  
  #################
  ### Add point ###
  #################
  while(blade.cost[length(blade.cost)] < 45.5){ # if total cost is less than the budget
    ### RNAmf ###
    predy <- predict(fit.RNAmf.next, X.test)$mu
    predsig2 <- predict(fit.RNAmf.next, X.test)$sig2
    
    ### RMSE ###
    blade.error <- c(blade.error, sqrt(mean((predy-y.test)^2))) # RMSE
    blade.crps <- c(blade.crps, mean(crps(y.test, predy, predsig2))) # CRPS
    if(alc.RNAmf$chosen$level == 1){
      blade.cost[length(blade.cost)+1] <- blade.cost[length(blade.cost)]+cost[1]
    }else{
      blade.cost[length(blade.cost)+1] <- blade.cost[length(blade.cost)]+sum(cost)
    }
    # print(blade.cost[length(blade.cost)])
    # print(blade.error[length(blade.error)])
    
    # if(blade.cost[length(blade.cost)] >= 45.5){break}
    
    ### update the next point ###
    alc.RNAmf <- ALC_RNAmf(Xref = X.test, Xcand = X.test, fit.RNAmf.next, cost = cost,
                           optim = TRUE, parallel = parallel, ncore = ncore)
    
    if (alc.RNAmf$chosen$level == 1) {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), alc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, blade1(alc.RNAmf$chosen$Xnext)),
                                        t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)],
                                        fit.RNAmf.next$fit2$y,
                                        kernel="sqex", constant=TRUE)
    } else {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), alc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, blade1(alc.RNAmf$chosen$Xnext)),
                                        rbind(t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)], alc.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit2$y, blade2(alc.RNAmf$chosen$Xnext)),
                                        kernel="sqex", constant=TRUE)
    }
    
    # save.image("Blade AL2.RData")
  }
  
  
  ### Save results ###
  costmatc2[[kk]] <- blade.cost
  rmsematc2[[kk]] <- blade.error
  crpsmatc2[[kk]] <- blade.crps
  
  time.each2[kk] <- proc.time()[3]- time.start
  # save.image("Blade ALC.RData")
}
costmatc2
rmsematc2
crpsmatc2
time.each2
