source("/Cokm;Internal_functions.R")
library(DiceKriging)
library(rgenoud)
library(mcmc)
library(MuFiCokriging)
library(lhs)
library(plgp)
library(matlabr)

one_step_cokm_varmax <- function(
    model,
    B,
    xpred,
    yreal = NULL,
    myfunctions, # no for blade
    niter,
    param.estim = FALSE,
    error.compute = FALSE, 
    error.LOO = TRUE, 
    ponderation = FALSE){
  
  Dnest <- model$Dnest
  y <- model$response
  
  nX <- dim(as.matrix(xpred))[1]
  nD <- Dnest$n
  d <- dim(as.matrix(xpred))[2]
  
  Dseq <- Dnest
  yseq <- y
  
  ypred <- predict(model, newdata=xpred, type="UK")
  res <- summary(model)
  
  rho <- res$Rho.Val[[1]]
  theta1 <- res$Cov.Val[[1]]
  theta2 <- res$Cov.Val[[2]]
  trend1 <- res$Trend.Val[[1]]
  trend2 <- res$Trend.Val[[2]]
  var1 <- res$Var.Val[[1]]
  var2 <- res$Var.Val[[2]]
  
  if(param.estim){
    coef.trend = NULL
    coef.cov = NULL
    coef.var = NULL
  }else{
    coef.trend = list(trend1,c(rho,trend2))
    coef.cov = list(theta1,theta2)
    coef.var = list(var1,var2)
  }
  
  if(error.compute){
    Q2 <- 1 - sum((yreal-ypred$mean)^2)/sum((yreal-mean(yreal))^2)
    RMSE <- mean((yreal-ypred$mean)^2)
    MaxAE <- max(abs(yreal-ypred$mean))	
  }
  if(error.LOO){
    Loo <- apply(matrix(1:length(y[[2]])),1,CVapply,model)
    Q2Loo <- 1 - sum(Loo^2)/sum((y[[2]]-mean(y[[2]]))^2)
    RMSELoo <- mean(Loo^2)
    MaxAELoo <- max(abs(Loo))		
    
    if(ponderation){
      #-- LOO variances
      varloo <- apply(matrix(1:length(y[[2]])),1,CVvarapply,model)
      #-- LOO errors
      e2loo <- (Loo)^2
      #-- Voronoi tessellation
      grid <- expand.grid(as.matrix(X)[,1],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,1])
      h <- (grid[,1]-grid[,2])^2
      if(d > 1){
        for(i in 2:d){
          grid <- expand.grid(as.matrix(X)[,i],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,i])
          h <-  h +(grid[,1]-grid[,2])^2
        }
      }
      matker <- matrix(h,nrow=nX, ncol= nD)
      indice <- max.col(-(matker))	
      #-- variance ponderation
      aj <- (e2loo/varloo)[indice]
      varaj <- adjustment(ypred$sig2,aj)
    }
  }
  
  ypredi <- ypred
  coutTot <- 0
  coutSave <- c()
  
  for(i in 1:niter){
    if(ponderation){
      ind.new <- which.max(varaj)
      xnew <- as.matrix(xpred)[ind.new,]
    }else{
      ind.new <- which.max(ypredi$sig2)
      xnew <- as.matrix(xpred)[ind.new,]
    }
    grid <- expand.grid(t(as.matrix(xnew))[,1],as.matrix(Dseq$PX)[,1])
    h <- (grid[,1]-grid[,2])^2
    if(d > 1){
      for(i in 2:d){
        grid <- expand.grid(t(as.matrix(xnew))[,i],as.matrix(Dseq$PX)[,i])
        h <-  h +(grid[,1]-grid[,2])^2
      }
    }
    exist <- FALSE
    if(min(h) < 1e-6){
      exist <- TRUE
      ind.exist <- which.min(h)
    }
    
    
    if(ponderation){
      #-- LOO variance
      varloo.delta <- apply(matrix(1:length(y[[2]])),1,CVvarapply1l,model)
      varloo.1 <- apply(matrix(1:length(y[[2]])),1,CVvarapplyLevel1,model)
      #-- LOO error
      Loo.delta <- apply(matrix(1:length(y[[2]])),1,CVapply1l,model)
      Loo.1 <- apply(matrix(1:length(y[[2]])),1,CVapplyLevel1,model)
      e2loo.delta <- (Loo.delta)^2
      e2loo.1 <- (Loo.1)^2
      #-- Voronoi tessellation 
      grid <- expand.grid(as.matrix(X)[,1],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,1])
      h <- (grid[,1]-grid[,2])^2
      if(d > 1){
        for(i in 2:d){
          grid <- expand.grid(as.matrix(X)[,i],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,i])
          h <-  h +(grid[,1]-grid[,2])^2
        }
      }
      matker <- matrix(h,nrow=nX, ncol= nD)
      indice <- max.col(-(matker))	
      #-- variance ponderation
      aj.delta <- (e2loo.delta/varloo.delta)[indice]
      varaj.delta <- adjustment(ypredi$sig2-rho^2*ypredi$varx[[1]],aj.delta)
      aj.1 <- (e2loo.1/varloo.1)[indice]
      varaj.1 <- adjustment(ypredi$varx[[1]],aj.1)
      
      ##-- criterion with adjusted variance
      IMSE.delta <- mean(varaj.delta)
      sig2.delta <- (varaj.delta)[ind.new]
      sig2.1 <- (varaj.1)[ind.new]
      
      if(sig2.delta < IMSE.delta){
        if(exist){"print: we reiterate at the same point"}else{
          Dnew <- rbind(Dseq$PX,xnew)
          row.names(Dnew) = NULL
          Dseq<- NestedDesign(Dnew, nlevel=2 , indices = Dseq$ind )
          
          # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
          
          xnew <- matrix(xnew, nrow=1)
          d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
          write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
          run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                            splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                            intern = TRUE)
          d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
          y1 <- c(y1, d2$V4)
          
          coutTot <- coutTot + 1
          coutSave <- c(coutSave,coutTot)
        }
      }else{
        critere <- prod(theta1)*rho^2*sig2.1/(prod(theta1)*rho^2*sig2.1+sig2.delta*prod(theta2))
        if(critere > 1/(1+B)){
          if(exist){"print: bug on reitere sur le meme point"}else{
            Dnew <- rbind(Dseq$PX,xnew)
            row.names(Dnew) = NULL
            Dseq<- NestedDesign(Dnew, nlevel=2 , indices = Dseq$ind )
            
            # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y1 <- c(y1, d2$V4)
            
            coutTot <- coutTot + 1
            coutSave <- c(coutSave,coutTot)
          }	
        }else{	
          if(exist){
            indice2 <- c(Dseq$ind[[1]],ind.exist)
            Dseq <- NestedDesign(Dseq$PX, nlevel=2 , indices = list(indice2) )
            
            # y2 <- c(y2,myfunctions[[2]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.0125) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y2 <- c(y2, d2$V4)
            
            coutTot <- coutTot + B + 1
            coutSave <- c(coutSave,coutTot)
          }else{
            indice2 <- c(Dseq$ind[[1]],(dim(Dseq$PX)[1]+1))
            Dnew <- rbind(Dseq$PX,xnew)
            row.names(Dnew) = NULL
            Dseq <- NestedDesign(Dnew, nlevel=2 , indices = list(indice2) )
            
            # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
            # y2 <- c(y2,myfunctions[[2]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y1 <- c(y1, d2$V4)
            d1 <- data.frame(xnew*0.5+0.25, 0.0125) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y2 <- c(y2, d2$V4)
            
            coutTot <- coutTot + B + 1
            coutSave <- c(coutSave,coutTot)
          }
        }
      }
    }else{
      IMSE.delta <- mean(ypredi$sig2-rho^2*ypredi$varx[[1]])
      sig2.delta <- (ypredi$sig2-rho^2*ypredi$varx[[1]])[ind.new]
      sig2.1 <- (ypredi$varx[[1]])[ind.new]
      
      if(sig2.delta < IMSE.delta){
        if(exist){"print: bug on reitere sur le meme point"}else{
          Dnew <- rbind(Dseq$PX,xnew)
          row.names(Dnew) = NULL
          Dseq<- NestedDesign(Dnew, nlevel=2 , indices = Dseq$ind )
          
          # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
          
          xnew <- matrix(xnew, nrow=1)
          d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
          write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
          run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                            splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                            intern = TRUE)
          d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
          y1 <- c(y1, d2$V4)
          
          coutTot <- coutTot + 1
          coutSave <- c(coutSave,coutTot)
        }
      }else{
        critere <- prod(theta1)*rho^2*sig2.1/(prod(theta1)*rho^2*sig2.1+sig2.delta*prod(theta2))
        if(critere > 1/(1+B)){
          if(exist){"print: bug on reitere sur le meme point"}else{
            Dnew <- rbind(Dseq$PX,xnew)
            row.names(Dnew) = NULL
            Dseq<- NestedDesign(Dnew, nlevel=2 , indices = Dseq$ind )
            
            # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y1 <- c(y1, d2$V4)
            
            coutTot <- coutTot + 1
            coutSave <- c(coutSave,coutTot)
          }	
        }else{
          if(exist){
            indice2 <- c(Dseq$ind[[1]],ind.exist)
            Dseq <- NestedDesign(Dseq$PX, nlevel=2 , indices = list(indice2) )
            
            # y2 <- c(y2,myfunctions[[2]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.0125) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y2 <- c(y2, d2$V4)
            
            coutTot <- coutTot + B + 1
            coutSave <- c(coutSave,coutTot)
          }else{
            indice2 <- c(Dseq$ind[[1]],(dim(Dseq$PX)[1]+1))
            Dnew <- rbind(Dseq$PX,xnew)
            row.names(Dnew) = NULL
            Dseq <- NestedDesign(Dnew, nlevel=2 , indices = list(indice2) )
            
            # y1 <- c(y1,myfunctions[[1]](t(as.matrix(xnew))))
            # y2 <- c(y2,myfunctions[[2]](t(as.matrix(xnew))))
            
            xnew <- matrix(xnew, nrow=1)
            d1 <- data.frame(xnew*0.5+0.25, 0.05) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y1 <- c(y1, d2$V4)
            d1 <- data.frame(xnew*0.5+0.25, 0.0125) # scale X to [-1,1]
            write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
            run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                              splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                              intern = TRUE)
            d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
            y2 <- c(y2, d2$V4)
            
            coutTot <- coutTot + B + 1
            coutSave <- c(coutSave,coutTot)
          }
        }
      }
    }
    
    modeli <- MuFicokm(
      formula = list(~1,~1), 
      MuFidesign = Dseq, 
      response = list(y1,y2), 
      nlevel = 2,
      covtype = "gauss",   #"matern5_2",
      coef.trend = coef.trend, 
      coef.cov = coef.cov, 
      coef.var = coef.var,
      nugget = 1e-06
    )
    
    ypredi <- predict(modeli, newdata=xpred, type="UK")
    if(error.compute){
      Q2 <- c(Q2,1 - sum((yreal-ypredi$mean)^2)/sum((yreal-mean(yreal))^2))
      RMSE <- c(RMSE,mean((yreal-ypredi$mean)^2))
      MaxAE <- c(MaxAE,max(abs(yreal-ypredi$mean)))
    }
    if(error.LOO){
      Loo <- apply(matrix(1:length(y[[2]])),1,CVapply,modeli)
      Q2Loo <- c(Q2Loo, 1 - sum(Loo^2)/sum((y[[2]]-mean(y[[2]]))^2))
      RMSELoo <- c(RMSELoo, mean(Loo^2))
      MaxAELoo <- c(MaxAELoo, max(abs(Loo)))		
      
      if(ponderation){
        #-- variance LOO
        varloo <- apply(matrix(1:length(y[[2]])),1,CVvarapply,modeli)
        #-- erreur LOO
        e2loo <- (Loo)^2
        #-- construction cellules voronoi  
        grid <- expand.grid(as.matrix(X)[,1],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,1])
        h <- (grid[,1]-grid[,2])^2
        if(d > 1){
          for(i in 2:d){
            grid <- expand.grid(as.matrix(X)[,i],as.matrix(Dseq$PX[Dseq$ind[[1]],])[,i])
            h <-  h +(grid[,1]-grid[,2])^2
          }
        }
        matker <- matrix(h,nrow=nX, ncol= nD)
        indice <- max.col(-(matker))	
        #-- ponderation
        aj <- (e2loo/varloo)[indice]
        varaj <- adjustment(ypredi$sig2,aj)
      }
    }
    if(param.estim){
      coef.trend = NULL
      coef.cov = NULL
      coef.var = NULL
    }else{
      coef.trend = list(trend1,c(rho,trend2))
      coef.cov = list(theta1,theta2)
      coef.var = list(var1,var2)
    }
  }
  if(error.compute){
    ValExt <- list(
      Q2 = Q2,
      RMSE = RMSE,
      MaxAE = MaxAE
    )
  }else{ValExt <- NULL}
  if(error.LOO){
    ValLoo <- list(
      Q2Loo = Q2Loo,
      RMSELoo = RMSELoo,
      MaxAELoo = MaxAELoo
    )
  }else{ValLoo <- NULL}
  
  
  return(list(Dseq = Dseq, 
              yseq = yseq, 
              modelseq = modeli,
              ypredseq = ypredi,
              D=Dnest,y=y,
              model = model, 
              ypred = ypred,
              ValExt = ValExt,
              ValLoo = ValLoo,
              CoutTot = coutTot,
              CoutSave = coutSave)
  )
}


eps <- sqrt(.Machine$double.eps)
crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

costmatco <- list(NA)
rmsematco <- list(NA)
crpsmatco <- list(NA)


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
  d1 <- data.frame(X1*0.5+0.25, rep(0.05, nrow(X1))) # scale X to [-1,1]
  write.csv(X1, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y1 <- d2$V4
  
  ### Y2 ###
  d1 <- data.frame(X2*0.5+0.25, rep(0.0125, nrow(X2))) # scale X to [-1,1]
  write.csv(X2, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y2 <- d2$V4
  
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
  print(sqrt(mean((predictions$mean-y.test)^2)))
  rmsematco[[kk]] <- sqrt(mean((predictions$mean-y.test)^2)) # Cokm
  crpsmatco[[kk]] <- mean(crps(y.test, predictions$mean, predictions$sig2)) # Cokm
  
  ## One point at-a-time sequential cokriging (see Section 3.1). ##
  B <- 6.85/2.25		#-B: ratio of computational costs between level 1 and 2
  
  for(i in 1:20){
    
    set.seed(kk)
    niter <- i
    
    cokm_varmax <- one_step_cokm_varmax(
      model = mymodel,
      B = B,
      xpred = X.test,
      yreal = y.test,
      # myfunctions = list(f1,f2),
      niter = niter,
      param.estim = TRUE,
      error.compute = FALSE,
      error.LOO = FALSE,
      ponderation = FALSE)
    
    ###
    print(sqrt(mean((cokm_varmax$ypredseq$mean-y.test)^2)))
    rmsematco[[kk]][i+1] <- sqrt(mean((cokm_varmax$ypredseq$mean-y.test)^2)) # Cokm
    crpsmatco[[kk]][i+1] <- mean(crps(y.test, cokm_varmax$ypredseq$mean, cokm_varmax$ypredseq$sig2)) # Cokm
    save.image("/Blade Cokm.RData")
  }
  
  costmatco[[kk]] <- 2.25*c(0,cokm_varmax$CoutSave)
  save.image("/Blade Cokm.RData")
  
}
costmatco
rmsematco
crpsmatco



