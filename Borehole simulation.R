### Borehole Example ###
borehole <- function(xx)
{
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

boreholelow <- function(xx)
{ 
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 5 * Tu * (Hu-Hl) #+ (Tu * Kw) # Tu * Kw is added
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1.5+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

output.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  borehole(x[1:8])
} 

outputlow.f <- function(x){
  factor_range <- list("rw" = c(0.05, 0.15), "r" = c(100, 50000),
                       "Tu" = c(63070, 115600), "Hu" = c(990, 1110),
                       "Tl" = c(63.1, 116), "Hl" = c(700, 820),
                       "L" = c(1120, 1680), "Kw" = c(9855, 12045))
  for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
  boreholelow(x[1:8])
} 

### training data ###
n1 <- 80; n2 <- 40; rep <- 100
d <- 8
### test data ###
n.test <- 1000

cat("Borehole function: n1 =", n1, ", n2 =", n2, ", rep =", rep, "\n")

result.borehole.rmse <- matrix(NA, rep, 3)
result.borehole.meancrps <- matrix(NA, rep, 3)
result.borehole.comptime <- matrix(NA, rep, 3)
colnames(result.borehole.rmse) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.borehole.meancrps) <- c("RNAmf", "Cokriging", "NARGP")
colnames(result.borehole.comptime) <- c("RNAmf", "Cokriging", "NARGP")


for(i in 1:rep) {
  set.seed(i)

  X1 <- maximinLHS(n1, d)
  X2 <- maximinLHS(n2, d)
  
  NestDesign <- NestedDesignBuild(design = list(X1,X2))
  
  X1 <- NestDesign$PX
  X2 <- ExtractNestDesign(NestDesign,2)
  
  y1 <- apply(X1,1,outputlow.f)
  y2 <- apply(X2,1,output.f)
  
  ### test data ###
  x <- maximinLHS(n.test, d)
  
  # Generate the filename using paste0
  filename <- paste0("RDSfile/file", i, ".rds") # change path
  saveRDS(list(X1=X1, X2=X2, Y1=y1, Y2=y2, Xtest=x, Ytest=apply(x,1,output.f)), file = filename)
  
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
                           lower=eps, #upper=0.1,
                           # coef.trend = list(0,c(0,0)),
                           response = list(y1,y2), nlevel = 2)
  pred.muficokm <- predict(fit.muficokm, x, "SK")
  toc.cokm <- proc.time()[3]
  
  result.borehole.rmse[i,1] <- sqrt(mean((predy-apply(x,1,output.f))^2)) 
  result.borehole.rmse[i,2] <- sqrt(mean((pred.muficokm$mean-apply(x,1,output.f))^2)) 
  # result.borehole.rmse[i,3] <- py$error
  
  result.borehole.meancrps[i,1] <- mean(crps(apply(x,1,output.f), predy, predsig2))
  result.borehole.meancrps[i,2] <- mean(crps(apply(x,1,output.f), pred.muficokm$mean, pred.muficokm$sig2)) 
  # result.borehole.meancrps[i,3] <- py$crps
  
  result.borehole.comptime[i,1] <- toc.RNAmf - tic.RNAmf
  result.borehole.comptime[i,2] <- toc.cokm - tic.cokm
  # result.borehole.comptime[i,3] <- py$ctime
  
  # boxplot(result.borehole.rmse[1:i,, drop=FALSE])
}

