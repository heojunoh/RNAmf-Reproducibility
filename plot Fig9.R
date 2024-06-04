# CRPS
resultmatc <- resultmatc2 <- resultmatc3 <- resultmatc4 <- resultmatco <- resultmatk <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r1 <- r2 <- r3 <- r4 <- rco <- rk <- c()
  
  for(i in 1:(length(costmatc[[j]])-1)){
    r1 <- c(r1, rep(crpsmatc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  if(ALC){
    for(i in 1:(length(costmatc2[[j]])-1)){
      r2 <- c(r2, rep(crpsmatc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
    }
  }
  for(i in 1:(length(costmatc3[[j]])-1)){
    r3 <- c(r3, rep(crpsmatc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  for(i in 1:(length(costmatc4[[j]])-1)){
    r4 <- c(r4, rep(crpsmatc4[[j]][i], costmatc4[[j]][i+1] - costmatc4[[j]][i]))
  }
  for(i in 1:(length(costmatco[[j]])-1)){
    rco <- c(rco, rep(crpsmatco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  for(i in 1:(length(costmatk[[j]])-1)){
    rk <- c(rk, rep(crpsmatk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  resultmatc[,j] <- c(r1, crpsmatc[[j]][length(costmatc[[j]])])[1:44]
  if(ALC) resultmatc2[,j] <- c(r2, crpsmatc2[[j]][length(costmatc2[[j]])])[1:44]
  resultmatc3[,j] <- c(r3, crpsmatc3[[j]][length(costmatc3[[j]])])[1:44]
  resultmatc4[,j] <- c(r4, crpsmatc4[[j]][length(costmatc4[[j]])])[1:44]
  resultmatco[,j] <- c(rco, crpsmatco[[j]][length(costmatco[[j]])])[1:44]
  resultmatk[,j] <- c(rk, crpsmatk[[j]][length(costmatk[[j]])])[1:44]
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
  perdalresult.crps <- data.frame(x=rep(37:80,6), 
                                  Strategy=factor(c(rep("ALM",44), rep("ALC",44), rep("ALMC",44), rep("ALD",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                                  levels=c("ALD","ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                  Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                  Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                  Min = c(resultminc, resultminc2, resultminc3, resultminc4, resultminco, resultmink))
  
}else{
  perdalresult.crps <- data.frame(x=rep(37:80,5), 
                                  Strategy=factor(c(rep("ALM",44), rep("ALMC",44), rep("ALD",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                                  levels=c("ALD","ALM","ALMC","Cokriging-CV","MR-SUR")),
                                  Mean = c(resultmeanc, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                                  Max = c(resultmaxc, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                                  Min = c(resultminc, resultminc3, resultminc4, resultminco, resultmink))
  
}

plotalperd.crps <- ggplot(perdalresult.crps, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank") + theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        plot.margin = margin(t = 0, r = 10, b = 0, l = 30),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))+
  scale_x_continuous(expand=c(0,0), limits = c(37, 80))+
  labs(x="Costs", y = "CRPS") 


# RMSE
resultmatc <- resultmatc2 <- resultmatc3 <- resultmatc4 <- resultmatco <- resultmatk <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r1 <- r2 <- r3 <- r4 <- rco <- rk <- c()
  
  for(i in 1:(length(costmatc[[j]])-1)){
    r1 <- c(r1, rep(rmsematc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  if(ALC){
    for(i in 1:(length(costmatc2[[j]])-1)){
      r2 <- c(r2, rep(rmsematc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
    }
  }
  for(i in 1:(length(costmatc3[[j]])-1)){
    r3 <- c(r3, rep(rmsematc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  for(i in 1:(length(costmatc4[[j]])-1)){
    r4 <- c(r4, rep(rmsematc4[[j]][i], costmatc4[[j]][i+1] - costmatc4[[j]][i]))
  }
  for(i in 1:(length(costmatco[[j]])-1)){
    rco <- c(rco, rep(rmsematco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  for(i in 1:(length(costmatk[[j]])-1)){
    rk <- c(rk, rep(rmsematk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  resultmatc[,j] <- c(r1, rmsematc[[j]][length(costmatc[[j]])])[1:44]
  if(ALC) resultmatc2[,j] <- c(r2, rmsematc2[[j]][length(costmatc2[[j]])])[1:44]
  resultmatc3[,j] <- c(r3, rmsematc3[[j]][length(costmatc3[[j]])])[1:44]
  resultmatc4[,j] <- c(r4, rmsematc4[[j]][length(costmatc4[[j]])])[1:44]
  resultmatco[,j] <- c(rco, rmsematco[[j]][length(costmatco[[j]])])[1:44]
  resultmatk[,j] <- c(rk, rmsematk[[j]][length(costmatk[[j]])])[1:44]
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
  perdalresult <- data.frame(x=rep(37:80,6), 
                             Strategy=factor(c(rep("ALM",44), rep("ALC",44), rep("ALMC",44), rep("ALD",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                             levels=c("ALD","ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                             Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                             Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                             Min = c(resultminc, resultminc2, resultminc3, resultminc4, resultminco, resultmink))
  
}else{
  perdalresult <- data.frame(x=rep(37:80,5), 
                             Strategy=factor(c(rep("ALM",44), rep("ALMC",44), rep("ALD",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                             levels=c("ALD","ALM","ALMC","Cokriging-CV","MR-SUR")),
                             Mean = c(resultmeanc, resultmeanc3, resultmeanc4, resultmeanco, resultmeank),
                             Max = c(resultmaxc, resultmaxc3, resultmaxc4, resultmaxco, resultmaxk),
                             Min = c(resultminc, resultminc3, resultminc4, resultminco, resultmink))
}

plotalperd.rmse <- ggplot(perdalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank") + theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        plot.margin = margin(t = 0, r = 30, b = 0, l = 10),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))+
  scale_x_continuous(expand=c(0,0), limits = c(37, 80))+
  labs(x="Costs", y = "RMSE") 

plotalperd <- ggarrange(plotalperd.rmse, NULL, plotalperd.crps, ncol=3, nrow=1, widths = c(1, -0.1, 1), common.legend = TRUE, legend="bottom")
plotalperd