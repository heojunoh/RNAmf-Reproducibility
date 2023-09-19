crps <- function(x, mu, sig2){ # The smaller, the better (0 to infinity)
  if(any(sig2==0)) sig2[sig2==0] <- eps
  -sqrt(sig2)*(1/sqrt(pi)-2*dnorm((x-mu)/sqrt(sig2))-(x-mu)/sqrt(sig2)*(2*pnorm((x-mu)/sqrt(sig2))-1))
}

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

eps <- sqrt(.Machine$double.eps)

rep <- 10
result.blade.rmse <- matrix(NA, rep, 2)
colnames(result.blade.rmse) <- c("RNAmf", "Cokriging")
result.blade.meancrps <- matrix(NA, rep, 2)
colnames(result.blade.meancrps) <- c("RNAmf", "Cokriging")
result.blade.comptime <- matrix(NA, rep, 2)
colnames(result.blade.comptime) <- c("RNAmf", "Cokriging")

n1 <- 20; n2 <- 10

### Test data ###
n <- 100
set.seed(1)
X.test <- maximinLHS(n, d)
d1 <- data.frame(X.test*0.5+0.25, rep(0.0125, nrow(X.test))) # scale X to [-1,1]
write.csv(X.test, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
                  splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                  intern = TRUE)
d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
y.test <- d2$V4
y.test
write.csv(y.test, "/Rmatlab_files/generate_text/ytestmin.txt", row.names=F)

for(i in 1:rep) {
  
  # i <- 10 # 1-10
  set.seed(i)
  
  ### Generate Input ###
  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  ### Y1 ###
  d1 <- data.frame(X1*0.5+0.25, rep(0.025, nrow(X1))) # scale X to [-1,1]
  write.csv(X1, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y1 <- d2$V4
  y1
  py_run_file("/Blade simulation 1.py")
  
  ### Y2 ###
  d1 <- data.frame(X2*0.5+0.25, rep(0.0125, nrow(X2))) # scale X to [-1,1]
  write.csv(X2, "/Rmatlab_files/generate_text/temp_to_X.txt", row.names=F)
  write.csv(d1, "/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  d2 <- read.table("/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  y2 <- d2$V4
  y2
  
  
  ### closed ###
  tic.closed <- proc.time()[3]
  fit.closed <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predRNAmf(fit.closed, X.test)$mu
  predsig2 <- predRNAmf(fit.closed, X.test)$sig2
  toc.closed <- proc.time()[3]
  
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=1,
                           # coef.trend = list(0,c(0,0),c(0,0)), 
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, X.test, "SK")
  toc.cokm <- proc.time()[3]
  
  py_run_file("/Blade simulation 2.py")
  
  ### RMSE ###
  result.blade.rmse[i,1] <- sqrt(mean((predy-y.test)^2)) # RNAmf
  result.blade.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-y.test)^2)) # Cokriging
  
  result.blade.meancrps[i,1] <- mean(crps(y.test, predy, predsig2)) # RNAmf
  result.blade.meancrps[i,2] <- mean(crps(y.test, pred.muficokm$mean, pred.muficokm$sig2)) # Cokriging
  
  result.blade.comptime[i,1] <- toc.closed - tic.closed
  result.blade.comptime[i,2] <- toc.cokm - tic.cokm
  
}

result.blade.rmse <- cbind(result.blade.rmse, NARGP=unlist(py$l2error))
result.blade.meancrps <- cbind(result.blade.meancrps, NARGP=unlist(py$meancrps))
result.blade.comptime <- cbind(result.blade.comptime, NARGP=unlist(py$comptime))
