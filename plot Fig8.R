# Computation time (sec.)
perd.comptime <- c(result.perd.comptime[,1], result.perd.comptime[,2], result.perd.comptime[,3])
branin.comptime <- c(result.branin.comptime[,1], result.branin.comptime[,2], result.branin.comptime[,3])
park.comptime <- c(result.park.comptime[,1], result.park.comptime[,2], result.park.comptime[,3])
borehole.comptime <- c(result.borehole.comptime[,1], result.borehole.comptime[,2], result.borehole.comptime[,3]) 
currin.comptime <- c(result.currin.comptime[,1], result.currin.comptime[,2], result.currin.comptime[,3])
franke.comptime <- c(result.franke.comptime[,1], result.franke.comptime[,2], result.franke.comptime[,3])

df.comptime <- data.frame(comptime=c(perd.comptime, branin.comptime, park.comptime, borehole.comptime, currin.comptime, franke.comptime), 
                          model=rep(model,6), 
                          Function=c(rep("Perdikaris",300),rep("Branin",300),rep("Park",300),rep("Borehole",300),rep("Currin",300),rep("Franke",300)), 
                          xx=c(rep("Perdikaris-RNAmf",100),rep("Perdikaris-CoKriging",100),rep("Perdikaris-NARGP",100),
                               rep("Branin-RNAmf",100),rep("Branin-CoKriging",100),rep("Branin-NARGP",100),
                               rep("Park-RNAmf",100), rep("Park-CoKriging",100), rep("Park-NARGP",100),
                               rep("Borehole-RNAmf",100),rep("Borehole-CoKriging",100),rep("Borehole-NARGP",100),
                               rep("Currin-RNAmf",100),rep("Currin-CoKriging",100),rep("Currin-NARGP",100),
                               rep("Franke-RNAmf",100),rep("Franke-CoKriging",100),rep("Franke-NARGP",100)) )
df.comptime$Function <- factor(df.comptime$Function , levels=c("Perdikaris", "Branin", "Park", "Borehole", "Currin", "Franke"))

numerical.time <- ggplot(df.comptime, aes(x=Function, y=comptime, fill=model, color=model)) + 
  geom_boxplot(alpha=0.5) + theme_ipsum() + ggtitle("") +
  theme(axis.title.x = element_blank(), legend.position="bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = margin(t = 0, r = 10, b = -0, l = 10),
        axis.title.y = element_text(size=16, margin = margin(r = 10), hjust=0.5),
        strip.text.x = element_blank(), text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Computational time (sec.)") + theme(panel.spacing = unit(0, "lines")) +
  geom_vline(xintercept = c(0.4,6.6), linetype="dotted") + facet_grid(. ~ model)
