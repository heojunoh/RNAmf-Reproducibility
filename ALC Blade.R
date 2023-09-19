library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(doParallel)
library(foreach)
library(RNAmf)
library(matlabr)

eps <- sqrt(.Machine$double.eps)
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

costmatc2 <- list(NA)
rmsematc2 <- list(NA)
crpsmatc2 <- list(NA)
time.each2 <- rep(0,10)
cost <- c(2.25, 6.85)

### setting ###
n1 <- 20; n2 <- 10
d <- 2              # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
Lmax <- 5
k <- NULL           
n.max <- 300        # the maximum number of sample size
log.fg <- TRUE
l <- 5

### test data ###
d <- 2
n <- 100
set.seed(1)
X.test <- maximinLHS(n, d)
d1 <- data.frame(X.test*0.5+0.25, rep(0.0125, nrow(X.test))) # scale X to [-1,1]
write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
write.csv(X.test, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                  splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                  intern = TRUE)
d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
y.test <- d2$V4 # MM=0.3 takes 10 minutes #


for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  print(kk)
  
  ### Generate Input ###
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  ### Y1 ###
  d1 <- data.frame(X1*0.5+0.25, rep(0.05, nrow(X1))) # scale X to [0.25,0.75]
  write.csv(X1, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y1 <- d2$V4
  
  ### Y2 ###
  d1 <- data.frame(X2*0.5+0.25, rep(0.0125, nrow(X2))) # scale X to [0.25,0.75]
  write.csv(X2, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y2 <- d2$V4
  
  ### RNAmf ###
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.RNAmf, X.test)$mu
  predsig2 <- predRNAmf(fit.RNAmf, X.test)$sig2
  
  ### RMSE ###
  print(sqrt(mean((predy-y.test)^2))) 
  
  
  blade.cost <- 0
  blade.error <- sqrt(mean((predy-y.test)^2))
  blade.crps <- mean(crps(y.test, predy, predsig2))
  
  Iselect <- ALC_two_level(X.test, fit.RNAmf, 100, cost, list(NULL, NULL), parallel=TRUE, ncore=10)
  
  
  #################
  ### Add point ###
  #################
  while(blade.cost[length(blade.cost)] < 45.5){ # if total cost is less than the budget
    ### closed ###
    predy <- predRNAmf(Iselect$fit, X.test)$mu
    predsig2 <- predRNAmf(Iselect$fit, X.test)$sig2
    
    ### RMSE ###
    blade.error <- c(blade.error, sqrt(mean((predy-y.test)^2))) # RMSE
    blade.crps <- c(blade.crps, mean(crps(y.test, predy, predsig2))) # CRPS
    if(Iselect$chosen$level == 1){
      blade.cost[length(blade.cost)+1] <- blade.cost[length(blade.cost)]+cost[1]
    }else{
      blade.cost[length(blade.cost)+1] <- blade.cost[length(blade.cost)]+sum(cost)
    }
    print(blade.cost[length(blade.cost)])
    print(blade.error[length(blade.error)])
    
    if(blade.cost[length(blade.cost)] >= 45.5){break}
    
    ### update the next point ###
    Iselect <- ALC_two_level(X.test, Iselect$fit, 100, cost, list(NULL, NULL), parallel=TRUE, ncore=10)
    save.image("/Blade AL2.RData")
  }
  
  
  ### Save results ###
  costmatc2[[kk]] <- blade.cost
  rmsematc2[[kk]] <- blade.error
  crpsmatc2[[kk]] <- blade.crps
  
  time.each2[kk] <- proc.time()[3]- time.start
  save.image("/Blade AL2.RData")
}
costmatc2
rmsematc2
crpsmatc2
time.each2



