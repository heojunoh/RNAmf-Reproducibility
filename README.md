Active Learning for a Recursive Non-Additive Emulator for Multi-Fidelity
Computer Experiments (Reproducibility)
================
September 18, 2023

This instruction aims to reproduce the results in the paper “*Active
Learning for a Recursive Non-Additive Emulator for Multi-Fidelity
Computer Experiments*”.

The following results are reproduced in this file

- Section 5.1: Figures 4, 5, and 6
- Section 5.2: Figures 7, 8, 9 and 10
- Section 6: Figure 12, 13 and 14

Section 5.2 and section 6 are run with 10 cores. You can change the
number of cores when you run. Followings are running time for each
section; Section 5.1; ~ 12 hours, Section 5.2; ~ 48 hours, Section 6; ~
12 hours. \[WARNING\] Reproducing results of Section 5.2 can take a long
time to run (more than a day). We recommend reproducing Section 5.1 and
Section 6 first before reproducing Section 5.2. Package “MuFiCokriging”
was removed from the CRAN repository. You can download “MuFiCokriging”
by following codes.

``` r
library(devtools)
install_github("cran/MuFiCokriging")
```

##### Step 0.1: load functions and packages

Before running .py file, get working directory by getwd.

``` r
library(lhs)
library(laGP)
library(plgp)
library(MuFiCokriging)
library(doParallel)
library(foreach)
library(reticulate)
library(rgenoud)
library(mcmc)
library(matlabr)
library(randtoolbox)
library(RColorBrewer)

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(hrbrthemes)
library(tidyverse)
library(gtable)
library(grid)
library(ggbreak)
library(patchwork)
extrafont::loadfonts()

source("/crps.R")
source("/Cokm;Internal_functions.R")    # This is from Le Gratiet and Cannamela (2015)
source("/Cokm;one_step_cokm_varmax.R")  # This is from Le Gratiet and Cannamela (2015)
py_install("GPy")
py_install("rpy2")
py_install("pandas")
py_install("matplotlib")
py_install("scipy")
py_install("time")

install_github("heojunoh/RNAmf")        # Download RNAmf
library(RNAmf)
```

##### Step 0.2: setting

``` r
eps <- sqrt(.Machine$double.eps) #small nugget for numeric stability
```

## Section 5.1:

##### Section 5.1: Reproducing Figure 4

This is reproducing plots of RMSE for synthetic functions. Each script
will run the simulation for the corresponding synthetic function. Python
codes may take long.

``` r
# Run 6 synthetic simulations
source("Perdikaris simulation.r")
source("Branin simulation.r")
source("Park simulation.r")
source("Borehole simulation.r")
source("Currin simulation.r")
source("Franke simulation.r")

model <- c(rep("RNAmf",100),rep("CoKriging",100),rep("NARGP",100))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))

# RMSE
perd.rmse <- c(result.perd.rmse[,1], result.perd.rmse[,2], result.perd.rmse[,3])

plot.perd.rmse <- ggplot(data.frame(perd.rmse, model), aes(x=model, y=perd.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Perdikaris") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


branin.rmse <- c(result.branin.rmse[,1], result.branin.rmse[,2], result.branin.rmse[,3])

plot.branin.rmse <- ggplot(data.frame(branin.rmse, model), aes(x=model, y=branin.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Branin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


park.rmse <- c(result.park.rmse[,1], result.park.rmse[,2], result.park.rmse[,3])

plot.park.rmse <- ggplot(data.frame(park.rmse, model), aes(x=model, y=park.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Park") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


borehole.rmse <- c(result.borehole.rmse[,1], result.borehole.rmse[,2], result.borehole.rmse[,3])

plot.borehole.rmse <- ggplot(data.frame(borehole.rmse, model), aes(x=model, y=borehole.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Borehole") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


currin.rmse <- c(result.currin.rmse[,1], result.currin.rmse[,2], result.currin.rmse[,3])

plot.currin.rmse <- ggplot(data.frame(currin.rmse, model), aes(x=model, y=currin.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Currin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


franke.rmse <- c(result.Franke.rmse[,1], result.Franke.rmse[,2], result.Franke.rmse[,3])

plot.franke.rmse <- ggplot(data.frame(franke.rmse, model), aes(x=model, y=franke.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Franke") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")

ggarrange(plot.perd.rmse, plot.branin.rmse, plot.park.rmse, 
          plot.borehole.rmse, plot.currin.rmse, plot.franke.rmse, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom") 
```

    ## Loading required package: ggplot2

<img src="figure/Synthetic RMSE.png" style="display: block; margin: auto;" />

##### Section 5.1: Reproducing Figure 5

This is reproducing plots of CRPS for synthetic functions.

``` r
# CRPS
perd.crps <- c(result.perd.meancrps[,1], result.perd.meancrps[,2], result.perd.meancrps[,3])

plot.perd.crps <- ggplot(data.frame(perd.crps, model), aes(x=model, y=perd.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Perdikaris") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


branin.crps <- c(result.branin.meancrps[,1], result.branin.meancrps[,2], result.branin.meancrps[,3])

plot.branin.crps <- ggplot(data.frame(branin.crps, model), aes(x=model, y=branin.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Branin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


park.crps <- c(result.park.meancrps[,1], result.park.meancrps[,2], result.park.meancrps[,3])

plot.park.crps <- ggplot(data.frame(park.crps, model), aes(x=model, y=park.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Park") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


borehole.crps <- c(result.borehole.meancrps[,1], result.borehole.meancrps[,2], result.borehole.meancrps[,3])

plot.borehole.crps <- ggplot(data.frame(borehole.crps, model), aes(x=model, y=borehole.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Borehole") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


currin.crps <- c(result.currin.meancrps[,1], result.currin.meancrps[,2], result.currin.meancrps[,3])

plot.currin.crps <- ggplot(data.frame(currin.crps, model), aes(x=model, y=currin.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Currin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


franke.crps <- c(result.Franke.meancrps[,1], result.Franke.meancrps[,2], result.Franke.meancrps[,3])

plot.franke.crps <- ggplot(data.frame(franke.crps, model), aes(x=model, y=franke.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Franke") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16, family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")

ggarrange(plot.perd.crps, plot.branin.crps, plot.park.crps, 
          plot.borehole.crps, plot.currin.crps, plot.franke.crps, 
          ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
```

<img src="figure/Synthetic CRPS.png" style="display: block; margin: auto;" />

##### Section 5.1: Reproducing Figure 6

This is reproducing plots of Computation time (sec.) for synthetic
functions.

``` r
# Computation time (sec.)
perd.comptime <- c(result.perd.comptime[,1], result.perd.comptime[,2], result.perd.comptime[,3])
branin.comptime <- c(result.branin.comptime[,1], result.branin.comptime[,2], result.branin.comptime[,3])
park.comptime <- c(result.park.comptime[,1], result.park.comptime[,2], result.park.comptime[,3])
borehole.comptime <- c(result.borehole.comptime[,1], result.borehole.comptime[,2], result.borehole.comptime[,3]) 
currin.comptime <- c(result.currin.comptime[,1], result.currin.comptime[,2], result.currin.comptime[,3])
franke.comptime <- c(result.franke.comptime[,1], result.franke.comptime[,2], result.franke.comptime[,3])

df.comptime <- data.frame(comptime=c(perd.comptime, branin.comptime, park.comptime, 
                                     borehole.comptime, currin.comptime, franke.comptime), 
                          model=rep(model,6), 
                          Function=c(rep("Perdikaris",300),rep("Branin",300),rep("Park",300),
                                     rep("Borehole",300),rep("Currin",300),rep("Franke",300)), 
                          xx=c(rep("Perdikaris-RNAmf",100),rep("Perdikaris-CoKriging",100),rep("Perdikaris-NARGP",100),
                               rep("Branin-RNAmf",100),rep("Branin-CoKriging",100),rep("Branin-NARGP",100),
                               rep("Park-RNAmf",100), rep("Park-CoKriging",100), rep("Park-NARGP",100),
                               rep("Borehole-RNAmf",100),rep("Borehole-CoKriging",100),rep("Borehole-NARGP",100),
                               rep("Currin-RNAmf",100),rep("Currin-CoKriging",100),rep("Currin-NARGP",100),
                               rep("Franke-RNAmf",100),rep("Franke-CoKriging",100),rep("Franke-NARGP",100))
                          )
df.comptime$Function <- factor(df.comptime$Function , levels=c("Perdikaris", "Branin", "Park", "Borehole", "Currin", "Franke"))


ggplot(df.comptime, 
       aes(x=Function, y=comptime, fill=model)) + 
  geom_boxplot(alpha=0.5)  +
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y = element_text(size=16, margin = margin(r = 10), hjust=0.5),
        strip.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "Computational time (sec.)")+
  theme(legend.position="none",
        panel.spacing = unit(0, "lines")) +
  geom_vline(xintercept = c(0.4,6.6), linetype="dotted")+
  facet_grid(. ~ model)
```

<img src="figure/Synthetic comptime.png" style="display: block; margin: auto;" />

## Section 5.2

##### Section 5.2: Reproducing Figure 7

This is reproducing plots of RMSE and CRPS for active learning on
Perdikaris function.

``` r
# Run 5 active learning strategies for Perdikaris function
source("ALM perdikaris.r")
source("ALC perdikaris.r")
source("ALMC perdikaris.r")
source("AL Cokm perdikaris.r")
source("AL MRSUR perdikaris.r")

# RMSE
resultmatc <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(rmsematc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, rmsematc[[j]][length(costmatc[[j]])])
  resultmatc[,j] <- r[1:44]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(rmsematc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, rmsematc2[[j]][length(costmatc2[[j]])])
  resultmatc2[,j] <- r[1:44]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r <- c(r, rep(rmsematc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  r <- c(r, rmsematc3[[j]][length(costmatc3[[j]])])
  resultmatc3[,j] <- r[1:44]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(rmsematco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, rmsematco[[j]][length(costmatco[[j]])])
  resultmatco[,j] <- r[1:44]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(rmsematk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, rmsematk[[j]][length(costmatk[[j]])])
  resultmatk[,j] <- r[1:44]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


perdalresult <- data.frame(x=rep(37:80,5), 
                                Strategy=factor(c(rep("ALM",44), rep("ALC",44), rep("ALMC",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                                levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                                Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                                Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalperd.rmse <- ggplot(perdalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_x_continuous(expand=c(0,0), limits = c(37, 80))+
  labs(x="Costs", y = "RMSE") 


# CRPS
resultmatc <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(crpsmatc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, crpsmatc[[j]][length(costmatc[[j]])])
  resultmatc[,j] <- r[1:44]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(crpsmatc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, crpsmatc2[[j]][length(costmatc2[[j]])])
  resultmatc2[,j] <- r[1:44]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r <- c(r, rep(crpsmatc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  r <- c(r, crpsmatc3[[j]][length(costmatc3[[j]])])
  resultmatc3[,j] <- r[1:44]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(crpsmatco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, crpsmatco[[j]][length(costmatco[[j]])])
  resultmatco[,j] <- r[1:44]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=44, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(crpsmatk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, crpsmatk[[j]][length(costmatk[[j]])])
  resultmatk[,j] <- r[1:44]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


perdalresult.crps <- data.frame(x=rep(37:80,5), 
                                Strategy=factor(c(rep("ALM",44), rep("ALC",44), rep("ALMC",44), rep("Cokriging-CV",44), rep("MR-SUR",44)), 
                                                levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                                Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                                Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalperd.crps <- ggplot(perdalresult.crps, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_x_continuous(expand=c(0,0), limits = c(37, 80))+
  labs(x="Costs", y = "CRPS") 

ggarrange(plotalperd.rmse, NULL, plotalperd.crps, ncol=3, nrow=1, widths = c(1, -0.1, 1), common.legend = TRUE, legend="bottom")
```

<img src="figure/Perd AL.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure 8

This is reproducing plots of final RMSE and proportion of low-fidelity
data for active learning on Perdikaris function.

``` r
# Final result
Strategy <- c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10))
Strategy <- factor(Strategy, levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

pperdfinal <- ggplot(data.frame(RMSE=c(resultmatc[51,], resultmatc2[51,], resultmatc3[51,], resultmatco[51,], resultmatk[51,]), Strategy=Strategy), aes(x=Strategy, y=RMSE, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Final RMSE") + theme(legend.position="bottom")


prop <- matrix(NA, ncol=5, nrow=10)

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc[[j]])-1){
    if(isTRUE(costmatc[[j]][i+1] - costmatc[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,1] <- (length(costmatc[[j]])-1)/(length(costmatc[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc2[[j]])-1){
    if(isTRUE(costmatc2[[j]][i+1] - costmatc2[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,2] <- (length(costmatc2[[j]])-1)/(length(costmatc2[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc3[[j]])-1){
    if(isTRUE(costmatc3[[j]][i+1] - costmatc3[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,3] <- (length(costmatc3[[j]])-1)/(length(costmatc3[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatco[[j]])-1){
    if(isTRUE(costmatco[[j]][i+1] - costmatco[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,4] <- (length(costmatco[[j]])-1)/(length(costmatco[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatk[[j]])-1){
    if(isTRUE(costmatk[[j]][i+1] - costmatk[[j]][i] == 1)){denom <- denom + 1}
  }
  prop[j,5] <- denom/(length(costmatk[[j]])-1)
}


df.propperd <- data.frame(Prop=c(prop[,1], prop[,2], prop[,3], prop[,4], prop[,5]), 
                          Strategy=c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10)))

df.propperd$Strategy <- factor(df.propperd$Strategy, levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

ppropperd <- ggplot(df.propperd, aes(x=Strategy, y=Prop, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5, size=20),
    axis.text.x = element_blank(),
    text=element_text(size=16),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Proportion of low-fidelity data")+
  theme(legend.position="none", panel.spacing = unit(0, "lines"))

(pperdfinal + ppropperd + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
```

<img src="figure/Perd AL prop.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure 9

This is reproducing plots of RMSE and CRPS for active learning on Park
function.

``` r
# Run 5 active learning strategies for Perdikaris function
source("ALM park.r")
source("ALC park.r")
source("ALMC park.r")
source("AL Cokm park.r")
source("AL MRSUR park.r")

# RMSE
resultmatc <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(rmsematc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, rmsematc[[j]][length(costmatc[[j]])])
  resultmatc[,j] <- r[1:31]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(rmsematc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, rmsematc2[[j]][length(costmatc2[[j]])])
  resultmatc2[,j] <- r[1:31]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r <- c(r, rep(rmsematc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  r <- c(r, rmsematc3[[j]][length(costmatc3[[j]])])
  resultmatc3[,j] <- r[1:31]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(rmsematco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, rmsematco[[j]][length(costmatco[[j]])])
  resultmatco[,j] <- r[1:31]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(rmsematk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, rmsematk[[j]][length(costmatk[[j]])])
  resultmatk[,j] <- r[1:31]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


parkalresult <- data.frame(x=rep(100:130,5), 
                           Strategy=factor(c(rep("ALM",31), rep("ALC",31), rep("ALMC",31), rep("Cokriging-CV",31), rep("MR-SUR",31)), 
                                           levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                           Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                           Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                           Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalpark.rmse <- ggplot(parkalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Costs", y = "RMSE") 


# CRPS
resultmatc <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(crpsmatc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, crpsmatc[[j]][length(costmatc[[j]])])
  resultmatc[,j] <- r[1:31]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(crpsmatc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, crpsmatc2[[j]][length(costmatc2[[j]])])
  resultmatc2[,j] <- r[1:31]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r <- c(r, rep(crpsmatc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  r <- c(r, crpsmatc3[[j]][length(costmatc3[[j]])])
  resultmatc3[,j] <- r[1:31]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(crpsmatco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, crpsmatco[[j]][length(costmatco[[j]])])
  resultmatco[,j] <- r[1:31]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=31, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(crpsmatk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, crpsmatk[[j]][length(costmatk[[j]])])
  resultmatk[,j] <- r[1:31]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


parkalresult.crps <- data.frame(x=rep(100:130,5), 
                           Strategy=factor(c(rep("ALM",31), rep("ALC",31), rep("ALMC",31), rep("Cokriging-CV",31), rep("MR-SUR",31)), 
                                           levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                           Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                           Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                           Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalpark.crps <- ggplot(parkalresult.crps, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Costs", y = "CRPS") 

ggarrange(plotalpark, NULL, plotalpark.crps, ncol=3, nrow=1, widths = c(1, -0.1, 1), common.legend = TRUE, legend="bottom")
```

<img src="figure/Park AL.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure 10

This is reproducing plots of final RMSE and proportion of low-fidelity
data for active learning on Park function.

``` r
# Final result
Strategy <- c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10))
Strategy <- factor(Strategy, levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

pparkfinal <- ggplot(data.frame(RMSE=c(resultmatc[31,], resultmatc2[31,], resultmatc3[31,], resultmatco[31,], resultmatk[31,]), Strategy=Strategy), aes(x=Strategy, y=RMSE, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Final RMSE") + theme(legend.position="bottom")


prop <- matrix(NA, ncol=5, nrow=10)

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc[[j]])-1){
    if(isTRUE(costmatc[[j]][i+1] - costmatc[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,1] <- (length(costmatc[[j]])-1)/(length(costmatc[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc2[[j]])-1){
    if(isTRUE(costmatc2[[j]][i+1] - costmatc2[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,2] <- (length(costmatc2[[j]])-1)/(length(costmatc2[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc3[[j]])-1){
    if(isTRUE(costmatc3[[j]][i+1] - costmatc3[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,3] <- (length(costmatc3[[j]])-1)/(length(costmatc3[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatco[[j]])-1){
    if(isTRUE(costmatco[[j]][i+1] - costmatco[[j]][i] == 4)){denom <- denom + 1}
  }
  prop[j,4] <- (length(costmatco[[j]])-1)/(length(costmatco[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatk[[j]])-1){
    if(isTRUE(costmatk[[j]][i+1] - costmatk[[j]][i] == 1)){denom <- denom + 1}
  }
  prop[j,5] <- denom/(length(costmatk[[j]])-1)
}


df.proppark <- data.frame(Prop=c(prop[,1], prop[,2], prop[,3], prop[,4], prop[,5]), 
                          Strategy=c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10)))


df.proppark$Strategy <- factor(df.proppark$Strategy , levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

pproppark <- ggplot(df.proppark, aes(x=Strategy, y=Prop, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5) + theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5, size=20),
    axis.text.x = element_blank(),
    text=element_text(size=16),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Proportion of low-fidelity data")+
  theme(legend.position="none", panel.spacing = unit(0, "lines"))

(pparkfinal + pproppark + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
```

<img src="figure/Park AL prop.png" style="display: block; margin: auto;" />

## Section 6

##### Section 6: Reproducing Figure 12

This is reproducing plots of RMSE and CRPS for Thermal Stress Analysis
of Jet Engine Turbine Blade data. This includes matlab and python codes.
Before you run, you should set the path on R, matlab, and python script.

``` r
# Run Blade simulations
source("Blade simulation.r")

model <- c(rep("RNAmf",10),rep("CoKriging",10),rep("NARGP",10))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))

blade.rmse <- c(result.blade.rmse[,1], result.blade.rmse[,2], result.blade.rmse[,3])

plot.blade.rmse <- ggplot(data.frame(blade.rmse, model), aes(x=model, y=blade.rmse, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "RMSE")+
  theme(legend.position="none")


blade.crps <- c(result.blade.meancrps[,1], result.blade.meancrps[,2], result.blade.meancrps[,3])

plot.blade.crps <- ggplot(data.frame(blade.crps, model), aes(x=model, y=blade.crps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "CRPS")+
  theme(legend.position="none")


blade.comptime <- c(result.blade.comptime[,1], result.blade.comptime[,2], result.blade.comptime[,3])

plot.blade.comptime <- ggplot(data.frame(blade.comptime, model), aes(x=model, y=blade.comptime, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "Computation time (sec.)")+
  theme(legend.position="none")


plot.blade <- ggarrange(plot.blade.rmse, NULL, plot.blade.crps, NULL, plot.blade.comptime, ncol=5, nrow=1, widths=c(1,-0.1,1,-0.1,1), common.legend = TRUE, legend="bottom")
```

<img src="figure/Blade comparison.png" style="display: block; margin: auto;" />

##### Section 6: Reproducing Figure 13

This is reproducing plots of RMSE and CRPS for active learning on Blade
data. Before you run, you should set the path on R and matlab script.

``` r
# Run 5 active learning strategies for Blade data
source("ALM blade.r")
source("ALC blade.r")
source("ALMC blade.r")
source("Cokm blade.r")
source("MRSUR blade.r")

# RMSE
resultmatc <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r[(costmatc[[j]][i]/0.05 + 1):(costmatc[[j]][i+1]/0.05)] <- rmsematc[[j]][i]
  }
  r <- c(r, rmsematc[[j]][length(rmsematc[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r[(costmatc2[[j]][i]/0.05 + 1):(costmatc2[[j]][i+1]/0.05)] <- rmsematc2[[j]][i]
  }
  r <- c(r, rmsematc2[[j]][length(rmsematc2[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc2[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r[(costmatc3[[j]][i]/0.05 + 1):(costmatc3[[j]][i+1]/0.05)] <- rmsematc3[[j]][i]
  }
  r <- c(r, rmsematc3[[j]][length(rmsematc3[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc3[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r[(costmatco[[j]][i]/0.05 + 1):(costmatco[[j]][i+1]/0.05)] <- rmsematco[[j]][i]
  }
  r <- c(r, rmsematco[[j]][length(rmsematco[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatco[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r[(costmatk[[j]][i]/0.05 + 1):(costmatk[[j]][i+1]/0.05)] <- rmsematk[[j]][i]
  }
  r <- c(r, rmsematk[[j]][length(rmsematk[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatk[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


bladealresult.rmse <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,5), 
                            Strategy=factor(c(rep("ALM",911), rep("ALC",911), rep("ALMC",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                            levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                            Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                            Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                            Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalblade.rmse <- ggplot(bladealresult.rmse, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_x_continuous(expand = c(0, 0), limits = c(113.5, 159))+
  labs(x="Costs", y = "RMSE") 


# CRPS
resultmatc <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r[(costmatc[[j]][i]/0.05 + 1):(costmatc[[j]][i+1]/0.05)] <- crpsmatc[[j]][i]
  }
  r <- c(r, crpsmatc[[j]][length(crpsmatc[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)


resultmatc2 <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r[(costmatc2[[j]][i]/0.05 + 1):(costmatc2[[j]][i+1]/0.05)] <- crpsmatc2[[j]][i]
  }
  r <- c(r, crpsmatc2[[j]][length(crpsmatc2[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc2[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)


resultmatc3 <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r[(costmatc3[[j]][i]/0.05 + 1):(costmatc3[[j]][i+1]/0.05)] <- crpsmatc3[[j]][i]
  }
  r <- c(r, crpsmatc3[[j]][length(crpsmatc3[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatc3[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)


resultmatco <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r[(costmatco[[j]][i]/0.05 + 1):(costmatco[[j]][i+1]/0.05)] <- crpsmatco[[j]][i]
  }
  r <- c(r, crpsmatco[[j]][length(crpsmatco[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatco[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)


resultmatk <- matrix(, nrow=length(seq(0, 45.5, 0.05)), ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r[(costmatk[[j]][i]/0.05 + 1):(costmatk[[j]][i+1]/0.05)] <- crpsmatk[[j]][i]
  }
  r <- c(r, crpsmatk[[j]][length(crpsmatk[[j]])])
  whichisna <- which(is.na(r))
  if(length(whichisna) > 0){
    for(k in 1:length(whichisna)){
      r[whichisna[k]] <- r[whichisna[k] - 1]
    }
  }
  if(length(r)==910) r[911] <- r[910]
  resultmatk[,j] <- r[1:length(seq(0, 45.5, 0.05))]
}

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


bladealresult.crps <- data.frame(x=rep(seq(0, 45.5, 0.05)+113.5,5), 
                                 Strategy=factor(c(rep("ALM",911), rep("ALC",911), rep("ALMC",911), rep("Cokriging-CV",911), rep("MR-SUR",911)), 
                                                 levels=c("ALM","ALC","ALMC","Cokriging-CV","MR-SUR")),
                                 Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                                 Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                                 Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalblade.crps <- ggplot(bladealresult.crps, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_x_continuous(expand = c(0, 0), limits = c(113.5, 159))+
  labs(x="Costs", y = "CRPS") 


ggarrange(plotalblade.rmse, NULL, plotalblade.crps, ncol=3, nrow=1, widths = c(1, -0.1, 1), common.legend = TRUE, legend="bottom")
```

<img src="figure/Blade AL.png" style="display: block; margin: auto;" />

##### Section 6: Reproducing Figure 14

This is reproducing plots of final RMSE and proportion of low-fidelity
data for active learning on Blade data.

``` r
# Final result
Strategy <- c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10))
Strategy <- factor(Strategy, levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

pbladefinal <- ggplot(data.frame(RMSE=c(resultmatc[911,],resultmatc2[911,],resultmatc3[911,],resultmatco[911,],resultmatk[911,]), Strategy=Strategy), aes(x=Strategy, y=RMSE, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Final RMSE")+
  theme(legend.position="bottom")


prop <- matrix(NA, ncol=5, nrow=10)

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc[[j]])-1){
    if(isTRUE(abs(costmatc[[j]][i+1] - costmatc[[j]][i] - (2.25+6.85)) < 1e-3)){denom <- denom + 1}
  }
  prop[j,1] <- (length(costmatc[[j]])-1)/(length(costmatc[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc2[[j]])-1){
    if(isTRUE(abs(costmatc2[[j]][i+1] - costmatc2[[j]][i] - (2.25+6.85)) < 1e-3)){denom <- denom + 1}
  }
  prop[j,2] <- (length(costmatc2[[j]])-1)/(length(costmatc2[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatc3[[j]])-1){
    if(isTRUE(abs(costmatc3[[j]][i+1] - costmatc3[[j]][i] - (2.25+6.85)) < 1e-3)){denom <- denom + 1}
  }
  prop[j,3] <- (length(costmatc3[[j]])-1)/(length(costmatc3[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatco[[j]])-1){
    if(isTRUE(abs(costmatco[[j]][i+1] - costmatco[[j]][i] - (2.25+6.85)) < 1e-3)){denom <- denom + 1}
  }
  prop[j,4] <- (length(costmatco[[j]])-1)/(length(costmatco[[j]])-1+denom)
}

for(j in 1:10){
  denom <- 0
  for(i in 1: length(costmatk[[j]])-1){
    if(isTRUE(abs(costmatk[[j]][i+1] - costmatk[[j]][i] - 2.25) < 1e-3)){denom <- denom + 1}
  }
  prop[j,5] <- denom/(length(costmatk[[j]])-1)
}


df.propblade <- data.frame(Prop=c(prop[,1], prop[,2], prop[,3], prop[,4], prop[,5]), 
Strategy=c(rep("ALM",10),rep("ALC",10),rep("ALMC",10),rep("Cokriging-CV",10),rep("MR-SUR",10)))

df.propblade$Strategy <- factor(df.propblade$Strategy , levels=c("ALM", "ALC", "ALMC", "Cokriging-CV", "MR-SUR"))

ppropblade <- ggplot(df.propblade, 
                    aes(x=Strategy, y=Prop, fill=Strategy, color=Strategy)) + 
  geom_boxplot(alpha=0.5)  +
  theme_bw() +  ggtitle("") +
  theme(
    axis.title.x = element_text(margin = margin(t = 10), hjust=0.5, size=20),
    axis.text.x = element_blank(),
    text=element_text(size=16),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0.5, size=20))+
  labs(x="", y = "Proportion of low-fidelity data")+
  theme(legend.position="none", panel.spacing = unit(0, "lines"))

(pbladefinal + ppropblade + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
```

<img src="figure/Blade AL prop.png" style="display: block; margin: auto;" />
