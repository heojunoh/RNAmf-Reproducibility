# CRPS
perd.crps <- c(result.perd.meancrps[,1], result.perd.meancrps[,2], result.perd.meancrps[,3])

plot.perd.crps <- ggplot(data.frame(perd.crps, model), aes(x=model, y=perd.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Perdikaris") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


branin.crps <- c(result.branin.meancrps[,1], result.branin.meancrps[,2], result.branin.meancrps[,3])

plot.branin.crps <- ggplot(data.frame(branin.crps, model), aes(x=model, y=branin.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Branin") + labs(x="", y = "")+
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


park.crps <- c(result.park.meancrps[,1], result.park.meancrps[,2], result.park.meancrps[,3])

plot.park.crps <- ggplot(data.frame(park.crps, model), aes(x=model, y=park.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Park") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


borehole.crps <- c(result.borehole.meancrps[,1], result.borehole.meancrps[,2], result.borehole.meancrps[,3])

plot.borehole.crps <- ggplot(data.frame(borehole.crps, model), aes(x=model, y=borehole.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Borehole") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


currin.crps <- c(result.currin.meancrps[,1], result.currin.meancrps[,2], result.currin.meancrps[,3])

plot.currin.crps <- ggplot(data.frame(currin.crps, model), aes(x=model, y=currin.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Currin") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))


franke.crps <- c(result.franke.meancrps[,1], result.franke.meancrps[,2], result.franke.meancrps[,3])

plot.franke.crps <- ggplot(data.frame(franke.crps, model), aes(x=model, y=franke.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_bw() + ggtitle("Franke") + labs(x="", y = "") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(), legend.position="none",
        plot.margin = margin(t = 10, r = 10, b = -10, l = 10),
        axis.title.y = element_blank(), panel.border = element_blank(),
        text=element_text(size=16), plot.title = element_text(hjust = 0.5, size=20))

numerical.crps <- ggarrange(plot.perd.crps, plot.branin.crps, plot.park.crps, 
                            plot.borehole.crps, plot.currin.crps, plot.franke.crps, 
                            ncol=3, nrow=2, common.legend = TRUE, legend="bottom")