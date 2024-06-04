model <- c(rep("RNAmf",10),rep("CoKriging",10),rep("NARGP",10))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))

blade.rmse <- c(result.blade.rmse[,1], result.blade.rmse[,2], result.blade.rmse[,3])

plot.blade.rmse <- ggplot(data.frame(blade.rmse, model), aes(x=model, y=blade.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_ipsum() + 
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 10),
        text=element_text(size=16), legend.position="none",
        plot.title = element_blank())+
  labs(x="", y = "RMSE")


blade.crps <- c(result.blade.meancrps[,1], result.blade.meancrps[,2], result.blade.meancrps[,3])

plot.blade.crps <- ggplot(data.frame(blade.crps, model), aes(x=model, y=blade.crps, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_ipsum() + 
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 10),
        text=element_text(size=16), legend.position="none",
        plot.title = element_blank())+
  labs(x="", y = "CRPS")


blade.comptime <- c(result.blade.comptime[,1], result.blade.comptime[,2], result.blade.comptime[,3])

plot.blade.comptime <- ggplot(data.frame(blade.comptime, model), aes(x=model, y=blade.comptime, fill=model)) + 
  geom_boxplot(alpha=0.5) + theme_ipsum() + 
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 0, l = 10),
        text=element_text(size=16), legend.position="none",
        plot.title = element_blank())+
  labs(x="", y = "Computation time (sec.)")


plot.blade <- ggarrange(plot.blade.rmse, NULL, plot.blade.crps, NULL, plot.blade.comptime, ncol=5, nrow=1, widths=c(1,-0.1,1,-0.1,1), common.legend = TRUE, legend="bottom")
plot.blade