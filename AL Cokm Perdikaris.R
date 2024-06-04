costmatco <- list(NA)
rmsematco <- list(NA)
crpsmatco <- list(NA)

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

#--- Initial experimental design sets
for(kk in 1:10){
  set.seed(kk) 
  
  cat("Perdikaris - Cokm: ", kk, "/ 10\n")
  
  X1 <- maximinLHS(n1, 1)
  X2 <- maximinLHS(n2, 1)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- f1(X1)
  y2 <- f2(X2)
  
  x <- seq(0,1,length.out=n.test)
  
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
    newdata = x,
    type="SK")
  ###
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-f2(x))^2)) # RMSE
  crpsmatco[[kk]] <- mean(crps(f2(x), predictions$mean, predictions$sig2)) # CRPS

  ## One point at-a-time sequential cokriging. ##
  B <- 3		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:40){ 
    ###
    set.seed(kk) 
    niter <- i 
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = x,
      yreal = f2(x),
      myfunctions = list(f1,f2),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-f2(x))^2)) # RMSE
    crpsmatco[[kk]][i+1] <- mean(crps(f2(x), cokm_varmax$ypredseq$mean, cokm_varmax$ypredseq$sig2)) # CRPS
  }

  costmatco[[kk]] <- c(0,cokm_varmax$CoutSave)
}
costmatco
rmsematco
crpsmatco
