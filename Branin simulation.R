### branin function ###
branin <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  
  (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+ 10
}

braninm <- function(xx)
{ 
  x1 <- xx[1]
  x2 <- xx[2]
  
  10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 + (10-5/(4*pi))*cos((x1+2))+ 10) + 2*(x1-0.5) - 3*(3*x2-1) - 1
}

braninl <- function(xx)
{ x1 <- xx[1]
x2 <- xx[2]

10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 + (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) + 2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
}

output.branin <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  branin(x[1:2])
} 

output.braninl <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  braninl(x[1:2])
} 

output.braninm <- function(x){
  factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
  
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  braninm(x[1:2])
} 


### training data ###
n1 <- 20; n2 <- 15; n3 <- 10; rep <- 100
d <- 2
### test data ###
n.test <- 1000

cat("Branin function: n1 =", n1, ", n2 =", n2, ", n3 =", n3, ", rep =", rep, "\n")


result.branin.rmse <- matrix(NA, rep, 3)
result.branin.meancrps <- matrix(NA, rep, 3)
result.branin.comptime <- matrix(NA, rep, 3)
colnames(result.branin.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.branin.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.branin.comptime) <- c("RNAmf", "Cokriging", "NARGP")

for(i in 1:rep) {
  set.seed(i)
  
  X1 <- maximinLHS(n1, 2)
  X2 <- maximinLHS(n2, 2)
  X3 <- maximinLHS(n3, 2)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  X3 <- ExtractNestDesign(NestDesign,3)
  
  y1 <- apply(X1,1,output.braninl)
  y2 <- apply(X2,1,output.braninm)
  y3 <- apply(X3,1,output.branin)
  
  ### test data ###
  x <- maximinLHS(n.test, d)
  
  # Generate the filename using paste0
  filename <- paste0("RDSfile/file", i, ".rds") # change path
  saveRDS(list(X1=X1, X2=X2, X3=X3, Y1=y1, Y2=y2, Y3=y3, Xtest=x, Ytest=apply(x,1,output.branin)), file = filename)
  
  ### RNAmf ###
  tic.RNAmf <- proc.time()[3]
  fit.RNAmf <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
  pred.RNAmf <- predict(fit.RNAmf, x)
  predy <- pred.RNAmf$mu
  predsig2 <- pred.RNAmf$sig2
  toc.RNAmf <- proc.time()[3]
  
  ### Cokriging ###
  tic.cokm <- proc.time()[3]
  fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                           lower=eps, upper=0.5,
                           # coef.trend = list(0,c(0,0),c(0,0)),
                           response = list(y1,y2,y3), nlevel = 3)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  

  result.branin.rmse[i,1] <- sqrt(mean((predy-apply(x,1,output.branin))^2)) 
  result.branin.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.branin))^2)) 
  # result.branin.rmse[i,3] <- py$error
  
  result.branin.meancrps[i,1] <- mean(crps(apply(x,1,output.branin), predy, predsig2))
  result.branin.meancrps[i,2] <- mean(crps(apply(x,1,output.branin), pred.muficokm$mean, pred.muficokm$sig2)) 
  # result.branin.meancrps[i,3] <- py$crps
  
  result.branin.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.branin.comptime[i,2] <- toc.cokm - tic.cokm
  # result.branin.comptime[i,3] <- py$ctime
  
  # boxplot(result.branin.rmse[1:i,, drop=FALSE])
}

