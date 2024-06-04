# Final result
if(ALC){
  Strategy <- c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("ALD",10),rep("Cokriging-CV",10),rep("MR-SUR",10))
  Strategy <- factor(Strategy, levels=c("ALD", "ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))
  
  pparkfinal <- ggplot(data.frame(RMSE=c(resultmatc[31,],resultmatc2[31,],resultmatc3[31,],resultmatc4[31,],resultmatco[31,],resultmatk[31,]), Strategy=Strategy), 
                       aes(x=Strategy, y=RMSE, fill=Strategy)) + 
    geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
    theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
          plot.margin = margin(t = 0, r = 30, b = 0, l = 10),
          panel.border = element_blank(), text=element_text(size=16),
          plot.title = element_text(hjust = 0.5, size=20))+
    labs(x="", y = "Final RMSE") + theme(legend.position="bottom")
}else{
  Strategy <- c(rep("ALM",10),rep("ALMC",10),rep("ALD",10),rep("Cokriging-CV",10),rep("MR-SUR",10))
  Strategy <- factor(Strategy, levels=c("ALD", "ALM", "ALMC", "Cokriging-CV", "MR-SUR"))
  
  pparkfinal <- ggplot(data.frame(RMSE=c(resultmatc[31,],resultmatc3[31,],resultmatc4[31,],resultmatco[31,],resultmatk[31,]), Strategy=Strategy), 
                       aes(x=Strategy, y=RMSE, fill=Strategy)) + 
    geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
    theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
          plot.margin = margin(t = 0, r = 30, b = 0, l = 10),
          panel.border = element_blank(), text=element_text(size=16),
          plot.title = element_text(hjust = 0.5, size=20))+
    labs(x="", y = "Final RMSE") + theme(legend.position="bottom")
}

prop <- matrix(NA, ncol=6, nrow=10)

for(j in 1:10){
  denom1 <- denom2 <- denom3 <- denom4 <- denomco <- denomk <- 0
  for(i in 1: length(costmatc[[j]])-1){
    if(isTRUE(costmatc[[j]][i+1] - costmatc[[j]][i] == 4)){denom1 <- denom1 + 1}
  }
  if(ALC){
    for(i in 1: length(costmatc2[[j]])-1){
      if(isTRUE(costmatc2[[j]][i+1] - costmatc2[[j]][i] == 4)){denom2 <- denom2 + 1}
    }
  }
  for(i in 1: length(costmatc3[[j]])-1){
    if(isTRUE(costmatc3[[j]][i+1] - costmatc3[[j]][i] == 4)){denom3 <- denom3 + 1}
  }
  for(i in 1: length(costmatc4[[j]])-1){
    if(isTRUE(costmatc4[[j]][i+1] - costmatc4[[j]][i] == 4)){denom4 <- denom4 + 1}
  }
  for(i in 1: length(costmatco[[j]])-1){
    if(isTRUE(costmatco[[j]][i+1] - costmatco[[j]][i] == 4)){denomco <- denomco + 1}
  }
  for(i in 1: length(costmatk[[j]])-1){
    if(isTRUE(costmatk[[j]][i+1] - costmatk[[j]][i] == 1)){denomk <- denomk + 1}
  }
  prop[j,1] <- (length(costmatc[[j]])-1)/(length(costmatc[[j]])-1+denom1)
  if(ALC) prop[j,2] <- (length(costmatc2[[j]])-1)/(length(costmatc2[[j]])-1+denom2)
  prop[j,3] <- (length(costmatc3[[j]])-1)/(length(costmatc3[[j]])-1+denom3)
  prop[j,4] <- (length(costmatc4[[j]])-1)/(length(costmatc4[[j]])-1+denom4)
  prop[j,5] <- (length(costmatco[[j]])-1)/(length(costmatco[[j]])-1+denomco)
  prop[j,6] <- denomk/(length(costmatk[[j]])-1)
}

if(ALC){
  df.proppark <- data.frame(Prop=c(prop[,1], prop[,2], prop[,3], prop[,4], prop[,5], prop[,6]), 
                            Strategy=c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("ALD",10),rep("Cokriging-CV",10),rep("MR-SUR",10)))
  df.proppark$Strategy <- factor(df.proppark$Strategy , levels=c("ALD", "ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))
}else{
  df.proppark <- data.frame(Prop=c(prop[,1], prop[,3], prop[,4], prop[,5], prop[,6]), 
                            Strategy=c(rep("ALM",10),rep("ALMC",10),rep("ALD",10),rep("Cokriging-CV",10),rep("MR-SUR",10)))
  df.proppark$Strategy <- factor(df.proppark$Strategy , levels=c("ALD", "ALM", "ALMC", "Cokriging-CV", "MR-SUR"))
}

pproppark <- ggplot(df.proppark, aes(x=Strategy, y=Prop, fill=Strategy)) + 
  geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5, size=20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        text=element_text(size=16), legend.position="none",
        panel.border = element_blank(), plot.margin = margin(t = 0, r = 10, b = 0, l = 30),
        plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Proportion of low-fidelity data")

plotalpark2 <- ggarrange(pparkfinal, pproppark, legend = "none")
plotalpark2