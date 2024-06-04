costmatco <- list(NA)
rmsematco <- list(NA)
crpsmatco <- list(NA)
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
  cat("Blade - Cokm: ", kk, "/ 10\n")
  
  ### Generate Input ###
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  ### Y1 ###
  y1 <- blade1(X1)
  
  ### Y2 ###
  y2 <- blade2(X2)
  
  #--- Multi-fidelity co-kriging model building
  mymodel <- MuFicokm(
    formula = list(~1,~1), 
    MuFidesign = NestDesign, 
    response = list(y1,y2), 
    lower=eps, upper=0.1,
    nlevel = 2,
    covtype = "gauss"
  )
  #--- Sequential design
  predictions <- predict(
    object = mymodel, 
    newdata = X.test,
    type="SK")
  ###
  # print(sqrt(mean((predictions$mean-y.test)^2)))
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-y.test)^2)) # Cokm
  crpsmatco[[kk]] <- mean(crps(y.test, predictions$mean, predictions$sig2)) # Cokm
  
  ## One point at-a-time sequential cokriging (see Section 3.1). ##
  B <- cost[2]/cost[1]		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:20){
    
    set.seed(kk)
    niter <- i
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = X.test,
      yreal = y.test,
      myfunctions = list(blade1, blade2),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    # print(sqrt(mean((cokm_varmax$ypredseq$mean-y.test)^2)))
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-y.test)^2)) # Cokm
    crpsmatco[[kk]][i+1] <- mean(crps(y.test, cokm_varmax$ypredseq$mean, cokm_varmax$ypredseq$sig2)) # Cokm
    # save.image("Blade Cokm.RData")
  }
  
  costmatco[[kk]] <- 2.25*c(0,cokm_varmax$CoutSave)
  # save.image("Blade Cokm.RData")
  
}
costmatco
rmsematco
crpsmatco
