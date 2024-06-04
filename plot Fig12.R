# CRPS
resultmatc <- resultmatc2 <- resultmatc3 <- resultmatc4 <- resultmatco <- resultmatk <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r1 <- r2 <- r3 <- r4 <- rco <- rk <- c()
  
  for(i in 1:(length(costmatc[[j]])-1)){
    r1[(costmatc[[j]][i]/0.05 + 1):(costmatc[[j]][i+1]/0.05)] <- crpsmatc[[j]][i]
  }
  r1 <- c(r1, crpsmatc[[j]][length(crpsmatc[[j]])])
  whichisna1 <- which(is.na(r1))
  if(length(whichisna1) > 0){
    for(k in 1:length(whichisna1)){
      r1[whichisna1[k]] <- r1[whichisna1[k] - 1]
    }
  }
  if(length(r1)==910) r1[911] <- r1[910]
  resultmatc[,j] <- r1[1:length(seq(0, 45.5, 0.05))]
  
  if(ALC){
    for(i in 1:(length(costmatc2[[j]])-1)){
      r2[(costmatc2[[j]][i]/0.05 + 1):(costmatc2[[j]][i+1]/0.05)] <- crpsmatc2[[j]][i]
    }
    r2 <- c(r2, crpsmatc2[[j]][length(crpsmatc2[[j]])])
    whichisna2 <- which(is.na(r2))
    if(length(whichisna2) > 0){
      for(k in 1:length(whichisna2)){
        r2[whichisna2[k]] <- r2[whichisna2[k] - 1]
      }
    }
    if(length(r2)==910) r2[911] <- r2[910]
    resultmatc2[,j] <- r2[1:length(seq(0, 45.5, 0.05))]
  }
  
  for(i in 1:(length(costmatc3[[j]])-1)){
    r3[(costmatc3[[j]][i]/0.05 + 1):(costmatc3[[j]][i+1]/0.05)] <- crpsmatc3[[j]][i]
  }
  r3 <- c(r3, crpsmatc3[[j]][length(crpsmatc3[[j]])])
  whichisna3 <- which(is.na(r3))
  if(length(whichisna3) > 0){
    for(k in 1:length(whichisna3)){
      r3[whichisna3[k]] <- r3[whichisna3[k] - 1]
    }
  }
  if(length(r3)==910) r3[911] <- r3[910]
  resultmatc3[,j] <- r3[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatc4[[j]])-1)){
    r4[(costmatc4[[j]][i]/0.05 + 1):(costmatc4[[j]][i+1]/0.05)] <- crpsmatc4[[j]][i]
  }
  r4 <- c(r4, crpsmatc4[[j]][length(crpsmatc4[[j]])])
  whichisna4 <- which(is.na(r4))
  if(length(whichisna4) > 0){
    for(k in 1:length(whichisna4)){
      r4[whichisna4[k]] <- r4[whichisna4[k] - 1]
    }
  }
  if(length(r4)==910) r4[911] <- r4[910]
  resultmatc4[,j] <- r4[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatc4[[j]])-1)){
    r4[(costmatc4[[j]][i]/0.05 + 1):(costmatc4[[j]][i+1]/0.05)] <- crpsmatc4[[j]][i]
  }
  r4 <- c(r4, crpsmatc4[[j]][length(crpsmatc4[[j]])])
  whichisna4 <- which(is.na(r4))
  if(length(whichisna4) > 0){
    for(k in 1:length(whichisna4)){
      r4[whichisna4[k]] <- r4[whichisna4[k] - 1]
    }
  }
  if(length(r4)==910) r4[911] <- r4[910]
  resultmatc4[,j] <- r4[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatco[[j]])-1)){
    rco[(costmatco[[j]][i]/0.05 + 1):(costmatco[[j]][i+1]/0.05)] <- crpsmatco[[j]][i]
  }
  rco <- c(rco, crpsmatco[[j]][length(crpsmatco[[j]])])
  whichisnaco <- which(is.na(rco))
  if(length(whichisnaco) > 0){
    for(k in 1:length(whichisnaco)){
      rco[whichisnaco[k]] <- rco[whichisnaco[k] - 1]
    }
  }
  if(length(rco)==910) rco[911] <- rco[910]
  resultmatco[,j] <- rco[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatk[[j]])-1)){
    rk[(costmatk[[j]][i]/0.05 + 1):(costmatk[[j]][i+1]/0.05)] <- crpsmatk[[j]][i]
  }
  rk <- c(rk, crpsmatk[[j]][length(crpsmatk[[j]])])
  whichisnak <- which(is.na(rk))
  if(length(whichisnak) > 0){
    for(k in 1:length(whichisnak)){
      rk[whichisnak[k]] <- rk[whichisnak[k] - 1]
    }
  }
  if(length(rk)==910) rk[911] <- rk[910]
  resultmatk[,j] <- rk[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)
if(ALC) {
  resultmeanc2 <- apply(resultmatc2, 1, mean)
  resultmaxc2 <- apply(resultmatc2, 1, max)
  resultminc2 <- apply(resultmatc2, 1, min)
}
resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)
resultmeanc4 <- apply(resultmatc4, 1, mean)
resultmaxc4 <- apply(resultmatc4, 1, max)
resultminc4 <- apply(resultmatc4, 1, min)
resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)
resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)

if(ALC){
  bladealresult.crps <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,6), 
                                   Strategy=factor(c(rep("ALM",911), rep("ALC",911), rep("ALMC",911), rep("ALD",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                                   levels=c("ALD","ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                   Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                   Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                   Min = c(resultminc, resultminc2, resultminc3, resultminc4, resultminco, resultmink))
}else{
  bladealresult.crps <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,5), 
                                   Strategy=factor(c(rep("ALM",911), rep("ALMC",911), rep("ALD",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                                   levels=c("ALD","ALM","ALMC","Cokriging-CV","MR-SUR")),
                                   Mean = c(resultmeanc, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                   Max = c(resultmaxc, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                   Min = c(resultminc, resultminc3, resultminc4, resultminco, resultmink))
}

plotalblade.crps <- ggplot(bladealresult.crps, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank") + theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 30),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  scale_y_continuous(limits = c(3, 25))+
  scale_x_continuous(expand = c(0, 0), limits = c(113.5, 159))+
  labs(x="Costs", y = "CRPS") 


# RMSE
resultmatc <- resultmatc2 <- resultmatc3 <- resultmatc4 <- resultmatco <- resultmatk <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  for(i in 1:(length(costmatc[[j]])-1)){
    r1[(costmatc[[j]][i]/0.05 + 1):(costmatc[[j]][i+1]/0.05)] <- rmsematc[[j]][i]
  }
  r1 <- c(r1, rmsematc[[j]][length(rmsematc[[j]])])
  whichisna1 <- which(is.na(r1))
  if(length(whichisna1) > 0){
    for(k in 1:length(whichisna1)){
      r1[whichisna1[k]] <- r1[whichisna1[k] - 1]
    }
  }
  if(length(r1)==910) r1[911] <- r1[910]
  resultmatc[,j] <- r1[1:length(seq(0, 45.5, 0.05))]
  
  if(ALC){
    for(i in 1:(length(costmatc2[[j]])-1)){
      r2[(costmatc2[[j]][i]/0.05 + 1):(costmatc2[[j]][i+1]/0.05)] <- rmsematc2[[j]][i]
    }
    r2 <- c(r2, rmsematc2[[j]][length(rmsematc2[[j]])])
    whichisna2 <- which(is.na(r2))
    if(length(whichisna2) > 0){
      for(k in 1:length(whichisna2)){
        r2[whichisna2[k]] <- r2[whichisna2[k] - 1]
      }
    }
    if(length(r2)==910) r2[911] <- r2[910]
    resultmatc2[,j] <- r2[1:length(seq(0, 45.5, 0.05))]
  }
  
  for(i in 1:(length(costmatc3[[j]])-1)){
    r3[(costmatc3[[j]][i]/0.05 + 1):(costmatc3[[j]][i+1]/0.05)] <- rmsematc3[[j]][i]
  }
  r3 <- c(r3, rmsematc3[[j]][length(rmsematc3[[j]])])
  whichisna3 <- which(is.na(r3))
  if(length(whichisna3) > 0){
    for(k in 1:length(whichisna3)){
      r3[whichisna3[k]] <- r3[whichisna3[k] - 1]
    }
  }
  if(length(r3)==910) r3[911] <- r3[910]
  resultmatc3[,j] <- r3[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatc4[[j]])-1)){
    r4[(costmatc4[[j]][i]/0.05 + 1):(costmatc4[[j]][i+1]/0.05)] <- rmsematc4[[j]][i]
  }
  r4 <- c(r4, rmsematc4[[j]][length(rmsematc4[[j]])])
  whichisna4 <- which(is.na(r4))
  if(length(whichisna4) > 0){
    for(k in 1:length(whichisna4)){
      r4[whichisna4[k]] <- r4[whichisna4[k] - 1]
    }
  }
  if(length(r4)==910) r4[911] <- r4[910]
  resultmatc4[,j] <- r4[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatco[[j]])-1)){
    rco[(costmatco[[j]][i]/0.05 + 1):(costmatco[[j]][i+1]/0.05)] <- rmsematco[[j]][i]
  }
  rco <- c(rco, rmsematco[[j]][length(rmsematco[[j]])])
  whichisnaco <- which(is.na(rco))
  if(length(whichisnaco) > 0){
    for(k in 1:length(whichisnaco)){
      rco[whichisnaco[k]] <- rco[whichisnaco[k] - 1]
    }
  }
  if(length(rco)==910) rco[911] <- rco[910]
  resultmatco[,j] <- rco[1:length(seq(0, 45.5, 0.05))]
  
  for(i in 1:(length(costmatk[[j]])-1)){
    rk[(costmatk[[j]][i]/0.05 + 1):(costmatk[[j]][i+1]/0.05)] <- rmsematk[[j]][i]
  }
  rk <- c(rk, rmsematk[[j]][length(rmsematk[[j]])])
  whichisnak <- which(is.na(rk))
  if(length(whichisnak) > 0){
    for(k in 1:length(whichisnak)){
      rk[whichisnak[k]] <- rk[whichisnak[k] - 1]
    }
  }
  if(length(rk)==910) rk[911] <- rk[910]
  resultmatk[,j] <- rk[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)
if(ALC){
  resultmeanc2 <- apply(resultmatc2, 1, mean)
  resultmaxc2 <- apply(resultmatc2, 1, max)
  resultminc2 <- apply(resultmatc2, 1, min)
}
resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)
resultmeanc4 <- apply(resultmatc4, 1, mean)
resultmaxc4 <- apply(resultmatc4, 1, max)
resultminc4 <- apply(resultmatc4, 1, min)
resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)
resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)

if(ALC){
  bladealresult.rmse <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,6), 
                                   Strategy=factor(c(rep("ALM",911), rep("ALC",911), rep("ALMC",911), rep("ALD",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                                   levels=c("ALD","ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                   Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                   Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                   Min = c(resultminc, resultminc2, resultminc3, resultminc4, resultminco, resultmink))
}else{
  bladealresult.rmse <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,5), 
                                   Strategy=factor(c(rep("ALM",911), rep("ALMC",911), rep("ALD",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                                   levels=c("ALD","ALM","ALMC","Cokriging-CV","MR-SUR")),
                                   Mean = c(resultmeanc, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                   Max = c(resultmaxc, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                   Min = c(resultminc, resultminc3, resultminc4, resultminco, resultmink))
}

plotalblade.rmse <- ggplot(bladealresult.rmse, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank") + theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        plot.margin = margin(t = 0, r = 30, b = 0, l = 10),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  scale_y_continuous(limits = c(5, 40))+
  scale_x_continuous(expand = c(0, 0), limits = c(113.5, 159))+
  labs(x="Costs", y = "RMSE") 


plotalblade <- ggarrange(plotalblade.rmse, NULL, plotalblade.crps, ncol=3, nrow=1, widths = c(1, -0.1, 1), common.legend = TRUE, legend="bottom")
plotalblade
