model <- c(rep("RNAmf",100),rep("CoKriging",100),rep("NARGP",100))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))

# RMSE
perd.rmse <- c(result.perd.rmse[,1], result.perd.rmse[,2], result.perd.rmse[,3])

plot.perd.rmse <- ggplot(data.frame(perd.rmse, model), aes(x=model, y=perd.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Perdikaris") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


branin.rmse <- c(result.branin.rmse[,1], result.branin.rmse[,2], result.branin.rmse[,3])

plot.branin.rmse <- ggplot(data.frame(branin.rmse, model), aes(x=model, y=branin.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Branin") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


park.rmse <- c(result.park.rmse[,1], result.park.rmse[,2], result.park.rmse[,3])

plot.park.rmse <- ggplot(data.frame(park.rmse, model), aes(x=model, y=park.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Park") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


borehole.rmse <- c(result.borehole.rmse[,1], result.borehole.rmse[,2], result.borehole.rmse[,3])

plot.borehole.rmse <- ggplot(data.frame(borehole.rmse, model), aes(x=model, y=borehole.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Borehole") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


currin.rmse <- c(result.currin.rmse[,1], result.currin.rmse[,2], result.currin.rmse[,3])

plot.currin.rmse <- ggplot(data.frame(currin.rmse, model), aes(x=model, y=currin.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Currin") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


franke.rmse <- c(result.franke.rmse[,1], result.franke.rmse[,2], result.franke.rmse[,3])

plot.franke.rmse <- ggplot(data.frame(franke.rmse, model), aes(x=model, y=franke.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Franke") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none", 
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

numerical.rmse <- ggarrange(plot.perd.rmse, plot.branin.rmse, plot.park.rmse, 
                            plot.borehole.rmse, plot.currin.rmse, plot.franke.rmse, 
                            ncol=3, nrow=2, common.legend = TRUE, legend="bottom") 