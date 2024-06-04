costmatco <- list(NA)
rmsematco <- list(NA)
crpsmatco <- list(NA)

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

#--- Initial experimental design sets
for(kk in 1:10){
  set.seed(kk) 
  
  cat("Park - Cokm: ", kk, "/ 10\n")
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)
  
  x <- maximinLHS(n.test, d)
  
  #--- Multi-fidelity co-kriging model building
  mymodel <- MuFicokm(
    formula = list(~1,~1), 
    MuFidesign = NestDesign, 
    response = list(y1,y2), 
    lower=eps, #upper=0.1,
    nlevel = 2,
    covtype = "gauss"
  )
  #--- Sequential design
  predictions <- predict(
    object = mymodel, 
    newdata = x,
    type="SK")
  ###
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-apply(x,1,park91a))^2)) # Cokm
  crpsmatco[[kk]] <- mean(crps(apply(x,1,park91a), predictions$mean, predictions$sig2)) # Cokm
  
  ## One point at-a-time sequential cokriging. ##
  B <- 3		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:25){ 
    ###
    set.seed(kk) 
    niter <- i 
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = x,
      yreal = park91a(x),
      myfunctions = list(park91alc, park91a),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-apply(x,1,park91a))^2)) # RMSE
    crpsmatco[[kk]][i+1] <- mean(crps(apply(x,1,park91a), cokm_varmax$ypredseq$mean, cokm_varmax$ypredseq$sig2)) # CRPS
  }
  
  costmatco[[kk]] <- c(0,cokm_varmax$CoutSave)
}
costmatco
rmsematco
crpsmatco
