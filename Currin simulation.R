### Currin Example ###
curretal88exp <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- 1 - exp(-1/(2*x2))
  fact2 <- 2300*x1^3 + 1900*x1^2 + 2092*x1 + 60
  fact3 <- 100*x1^3 + 500*x1^2 + 4*x1 + 20
  
  y <- fact1 * fact2/fact3
  return(y)
}

curretal88explc <- function(xx)
{
  x1 <- xx[1]
  x2 <- xx[2]
  
  maxarg <- max(c(0, x2-1/20))
  
  yh1 <- curretal88exp(c(x1+1/20, x2+1/20))
  yh2 <- curretal88exp(c(x1+1/20, maxarg))
  yh3 <- curretal88exp(c(x1-1/20, x2+1/20))
  yh4 <- curretal88exp(c(x1-1/20, maxarg))
  
  y <- (yh1 + yh2 + yh3 + yh4) / 4
  return(y)
}

### training data ###
n1 <- 20; n2 <- 10; rep <- 100
d <- 2
### test data ###
n.test <- 1000

cat("Currin function: n1 =", n1, ", n2 =", n2, ", rep =", rep, "\n")


result.currin.rmse <- matrix(NA, rep, 3)
result.currin.meancrps <- matrix(NA, rep, 3)
result.currin.comptime <- matrix(NA, rep, 3)
colnames(result.currin.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.currin.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.currin.comptime) <- c("RNAmf", "Cokriging", "NARGP")


for(i in 1:rep) {
  set.seed(i)

  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,curretal88explc)
  y2 <- apply(X2,1,curretal88exp)
  
  ### test data ###
  x <- maximinLHS(n.test, d)
  
  # Generate the filename using paste0
  filename <- paste0("RDSfile/file", i, ".rds") # change path
  saveRDS(list(X1=X1, X2=X2, Y1=y1, Y2=y2, Xtest=x, Ytest=apply(x,1,curretal88exp)), file = filename)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy <- pred.RNAmf$mu
  predsig2 <- pred.RNAmf$sig2
  toc.RNAmf <- proc.time()[3]
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  
  result.currin.rmse[i,1] <- sqrt(mean((predy-apply(x,1,curretal88exp))^2)) 
  result.currin.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,curretal88exp))^2)) 
  # result.currin.rmse[i,3] <- py$error
  
  result.currin.meancrps[i,1] <- mean(crps(apply(x,1,curretal88exp), predy, predsig2))
  result.currin.meancrps[i,2] <- mean(crps(apply(x,1,curretal88exp), pred.muficokm$mean, pred.muficokm$sig2)) 
  # result.currin.meancrps[i,3] <- py$crps
  
  result.currin.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.currin.comptime[i,2] <- toc.cokm - tic.cokm
  # result.currin.comptime[i,3] <- py$ctime
  
  # boxplot(result.currin.rmse[1:i,, drop=FALSE])
}

