cost <- NULL
d <- 2              # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
Lmax <- 5
k <- NULL           # k: the parameter for the Matern kernel function (NULL means it'll be estimated by LOOCV) 
n.max <- 300        # the maximum number of sample size
log.fg <- TRUE
l <- 5

blade1 <- function(xx){
  d1 <- data.frame(xx*0.5+0.25, rep(0.05, nrow(xx))) # scale X to [-1,1]
  write.csv(xx, "Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y <- d2$V4
  y
}

blade2 <- function(xx){
  d1 <- data.frame(xx*0.5+0.25, rep(0.0125, nrow(xx))) # scale X to [-1,1]
  write.csv(xx, "Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y <- d2$V4
  y
}


rep <- 10
result.blade.rmse <- matrix(NA, rep, 3)
result.blade.meancrps <- matrix(NA, rep, 3)
result.blade.comptime <- matrix(NA, rep, 3)
colnames(result.blade.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.blade.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.blade.comptime) <- c("RNAmf", "Cokriging", "NARGP")

n1 <- 20; n2 <- 10

### Test data ###
n <- 100
set.seed(1)
X.test <- maximinLHS(n, d)
y.test <- blade2(X.test)

for(i in 1:rep) {
  
  set.seed(i)
  
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
  
  
  eps <- sqrt(.Machine$double.eps)
  
  saveRDS(list(X1=X1, X2=X2, Y1=y1, Y2=y2, Xtest=X.test, Ytest=y.test), file = "tmp_data_blade.rds")
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predRNAmf(fit.RNAmf, X.test)
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
  
  
  ### NARGP ###
  py <- py_run_file("python code/Blade.py")
  
  ### RMSE ###
  result.blade.rmse[i,1] <- sqrt(mean((predy-y.test)^2)) # RNAmf
  result.blade.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) # Cokriging
  result.blade.rmse[i,3] <- py$error
  
  result.blade.meancrps[i,1] <- mean(crps(y.test, predy, predsig2)) # RNAmf
  result.blade.meancrps[i,2] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  result.blade.meancrps[i,3] <- py$crps
  
  result.blade.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.blade.comptime[i,2] <- toc.cokm - tic.cokm
  result.blade.comptime[i,3] <- py$ctime
}

