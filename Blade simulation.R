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
n1 <- 20; n2 <- 10; rep <- 10
d <- 2
### test data ###
n.test <- 100

cat("Blade application: n1 =", n1, ", n2 =", n2, ", rep =", rep, "\n")

result.blade.rmse <- matrix(NA, rep, 3)
result.blade.meancrps <- matrix(NA, rep, 3)
result.blade.comptime <- matrix(NA, rep, 3)
colnames(result.blade.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.blade.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.blade.comptime) <- c("RNAmf", "Cokriging", "NARGP")


### Test data ###
set.seed(1)
X.test <- maximinLHS(n.test, d)
y.test <- blade2(X.test)

for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  ### Y1 ###
  y1 <- blade1(X1)
  
  ### Y2 ###
  y2 <- blade2(X2)
  
  # Generate the filename using paste0
  filename <- paste0("RDSfile/file", i, ".rds") # change path
  saveRDS(list(X1=X1, X2=X2, Y1=y1, Y2=y2, Xtest=X.test, Ytest=y.test), file = filename)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, X.test)
  predy <- pred.RNAmf$mu
  predsig2 <- pred.RNAmf$sig2
  toc.RNAmf <- proc.time()[3]
  
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=1,
                           # coef.trend = list(0,c(0,0),c(0,0)), 
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, X.test, "SK")
  toc.cokm <- proc.time()[3]
  
  
  ### RMSE ###
  result.blade.rmse[i,1] <- sqrt(mean((predy-y.test)^2)) # RNAmf
  result.blade.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) # Cokriging
  # result.blade.rmse[i,3] <- py$error # NARGP
  
  result.blade.meancrps[i,1] <- mean(crps(y.test, predy, predsig2)) # RNAmf
  result.blade.meancrps[i,2] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  # result.blade.meancrps[i,3] <- py$crps # NARGP
  
  result.blade.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.blade.comptime[i,2] <- toc.cokm - tic.cokm
  # result.blade.comptime[i,3] <- py$ctime # NARGP
}

