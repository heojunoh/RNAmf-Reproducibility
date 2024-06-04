costmatc4 <- list(NA)
rmsematc4 <- list(NA)
crpsmatc4 <- list(NA)
time.each <- rep(0,10)
cost <- c(1,3)
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

for(kk in 1:10){
  time.start <- proc.time()[3]
  set.seed(kk)
  
  cat("Park - ALD: ", kk, "/ 10\n")

  X <- NestedX(c(n1, n2), d)
  X1 <- X[[1]]
  X2 <- X[[2]]

  y1 <- apply(X1,1,park91alc)
  y2 <- apply(X2,1,park91a)

  ### test data ###
  x <- maximinLHS(n.test, d)


  ### RNAmf ###
  fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
  predy <- predict(fit.RNAmf, x)$mu
  predsig2 <- predict(fit.RNAmf, x)$sig2

  ### RMSE ###
  # print(sqrt(mean((predy-apply(x,1,park91a))^2))) # closed form

  park.cost <- 0
  park.error <- sqrt(mean((predy-apply(x,1,park91a))^2))
  park.crps <- mean(crps(apply(x,1,park91a), predy, predsig2))

  ald.RNAmf <- ALD_RNAmf(Xcand = x, fit.RNAmf, cost = cost, optim = TRUE, parallel = parallel, ncore=ncore)

  if (ald.RNAmf$chosen$level == 1) {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), ald.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, apply(ald.RNAmf$chosen$Xnext,1,park91alc)),
                                      t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)],
                                      fit.RNAmf$fit2$y,
                                      kernel="sqex", constant=TRUE)
  } else {
    fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf$fit1$X)*attr(fit.RNAmf$fit1$X,"scaled:scale")+attr(fit.RNAmf$fit1$X,"scaled:center")), ald.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit1$y, apply(ald.RNAmf$chosen$Xnext,1,park91alc)),
                                      rbind(t(t(fit.RNAmf$fit2$X)*attr(fit.RNAmf$fit2$X,"scaled:scale")+attr(fit.RNAmf$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)], ald.RNAmf$chosen$Xnext),
                                      c(fit.RNAmf$fit2$y, apply(ald.RNAmf$chosen$Xnext,1,park91a)),
                                      kernel="sqex", constant=TRUE)
  }

  #################
  ### Add point ###
  #################
  while(park.cost[length(park.cost)] < 30){ # if total cost is less than the budget
    ### RNAmf ###
    predy <- predict(fit.RNAmf.next, x)$mu
    predsig2 <- predict(fit.RNAmf.next, x)$sig2

    ### RMSE ###
    park.error <- c(park.error, sqrt(mean((predy-apply(x,1,park91a))^2))) # RMSE
    park.crps <- c(park.crps, mean(crps(apply(x,1,park91a), predy, predsig2))) # CRPS
    if(ald.RNAmf$chosen$level == 1){
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+cost[1]
    }else{
      park.cost[length(park.cost)+1] <- park.cost[length(park.cost)]+sum(cost)
    }
    # print(park.cost[length(park.cost)])
    # print(park.error[length(park.error)])

    # if(park.cost[length(park.cost)] >= 30){break}

    ### update the next point ###
    ald.RNAmf <- ALD_RNAmf(Xcand = x, fit.RNAmf.next, cost = cost, optim = TRUE, parallel = parallel, ncore = ncore)

    if (ald.RNAmf$chosen$level == 1) {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), ald.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, apply(ald.RNAmf$chosen$Xnext,1,park91alc)),
                                        t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)],
                                        fit.RNAmf.next$fit2$y,
                                        kernel="sqex", constant=TRUE)
    } else {
      fit.RNAmf.next <- RNAmf_two_level(rbind(t(t(fit.RNAmf.next$fit1$X)*attr(fit.RNAmf.next$fit1$X,"scaled:scale")+attr(fit.RNAmf.next$fit1$X,"scaled:center")), ald.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit1$y, apply(ald.RNAmf$chosen$Xnext,1,park91alc)),
                                        rbind(t(t(fit.RNAmf.next$fit2$X)*attr(fit.RNAmf.next$fit2$X,"scaled:scale")+attr(fit.RNAmf.next$fit2$X,"scaled:center"))[,-ncol(fit.RNAmf$fit2$X)], ald.RNAmf$chosen$Xnext),
                                        c(fit.RNAmf.next$fit2$y, apply(ald.RNAmf$chosen$Xnext,1,park91a)),
                                        kernel="sqex", constant=TRUE)
    }

    # save.image("/Park AL4 1,3.RData")
  }

  ### Save results ###
  costmatc4[[kk]] <- park.cost
  rmsematc4[[kk]] <- park.error
  crpsmatc4[[kk]] <- park.crps

  time.each[kk] <- proc.time()[3]- time.start
  # save.image("Park ALD.RData")
}
costmatc4
rmsematc4
crpsmatc4
time.each
