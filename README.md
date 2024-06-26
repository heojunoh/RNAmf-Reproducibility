Active Learning for a Recursive Non-Additive Emulator for Multi-Fidelity
Computer Experiments (Reproducibility)
================
Junoh Heo, Chih-Li Sung
Jun 3, 2024

This instruction aims to reproduce the results in the paper “*Active
Learning for a Recursive Non-Additive Emulator for Multi-Fidelity
Computer Experiments*”.

The following results are reproduced in this file

- Section 5.1: Figures 7, S14, and 8
- Section 5.2: Figures 9, 10, S15, and S16
- Section 6: Figure 11, 12 and 13

The approximate running times for each section are as follows:

- Section 5.1: ~9 hours
- Section 5.2: ~48 hours
- Section 6: ~12 hours

**$$WARNING$$** Reproducing the results of Section 5.2 can take a
considerable amount of time (more than a day). We recommend reproducing
the results of Section 5.1 and Section 6 first before attempting Section
5.2.

The package `MuFiCokriging` has been removed from the CRAN repository.
You can download `MuFiCokriging` using the following code:

``` r
library(devtools)
install_github("cran/MuFiCokriging")
```

##### Step 0.1: setting

As mentioned earlier, reproducing all the results can be time-consuming.
Therefore, we have included options in the settings to allow users to
skip the most resource-intensive parts. Specifically, the logical object
`NARGP` determines whether you want to run the NARGP method for
comparison, which is implemented in Python. This method involves
sophisticated Monte Carlo (MC) approximations and can take a long time
to run. Additionally, installing the required Python packages can be
complex and time-consuming. Therefore, we recommend setting
**NARGP=FALSE** to avoid running this method during your initial
attempt.

Similarly, the logical object `ALC` determines whether you want to run
the active learning strategy “ALC”, which also involves complex MC
approximations, taking a long time to run. We recommend setting
**ALC=FALSE** to avoid running this strategy during your initial
attempt.

Section 5.2 and Section 6 are executed using 10 cores for parallel
computing. You can adjust the number of cores when you run the code.

``` r
setwd("/RNAmf-Reproducibility-main/") # set your working directory
eps <- sqrt(.Machine$double.eps) # small nugget for numeric stability
NARGP <- TRUE # do you want to run NARGP method
ALC <- TRUE   # do you want to run ALC strategy
parallel <- TRUE # do you want to run parallel for Section 5.2
ncore <- 10 # how many cores?
```

##### Step 0.2: load functions and packages

Before running .py files, you may need to set a proper working
environment. *This document is written based on Apple silicon (M3) Mac.
Different settings may be required for Intel Macs, Windows, or other
operating systems. Please carefully check.*

``` r
library(lhs)
library(plgp)
library(MuFiCokriging)
library(doParallel)
library(matlabr)
library(ggplot2)
library(ggpubr)
library(hrbrthemes)

source("crps.R") # CRPS evaluation

# This is for running cokriging, originally from Le Gratiet and Cannamela (2015)
source("Cokm;Internal_functions.R")   
source("Cokm;one_step_cokm_varmax.R") 
# You may need to install packages using pip install on terminal or command. Sometimes it requires to restart R/Rstudio after packages installed.
if(NARGP){
  library(reticulate)
  py_install("GPy==1.10.0")
  py_install("numpy==1.24.0")
  py_install("pandas==2.2.1")
  py_install("scipy==1.13.0")
  py_install("rpy2==3.5.11")
}

install.packages("RNAmf_0.1.2.tar.gz", repos = NULL, type = "source") # install RNAmf package
library(RNAmf)
```

## Section 5.1:

##### Section 5.1: Reproducing Figure 7

This reproduces the RMSE results for synthetic functions. Each script
will run the simulation for the corresponding synthetic function. These
scripts save the training samples with names `file1.rds` to
`file100.rds` in the `RDSfile` folder for fair comparison with NARGP.
Make sure to run the code one synthetic example at a time to avoid
overwriting files. Note that you may encounter the error *Error in
py_run_file_impl(file, local, convert) : numpy.linalg.LinAlgError: SVD
did not converge* when you run `Franke.py` for NARGP. We return `NA` for
those cases.

``` r
# Run 6 synthetic simulations
if(NARGP){
  source("Perdikaris simulation.R")
  py <- py_run_file("python code/Perdikaris.py")
  result.perd.rmse[,3] <- unlist(py$error)
  result.perd.meancrps[,3] <- unlist(py$crps)
  result.perd.comptime[,3] <- unlist(py$ctime)
    
  source("Branin simulation.R")
  py <- py_run_file("python code/Branin.py")
  result.branin.rmse[,3] <- unlist(py$error)
  result.branin.meancrps[,3] <- unlist(py$crps)
  result.branin.comptime[,3] <- unlist(py$ctime)
    
  source("Park simulation.R")
  py <- py_run_file("python code/Park.py")
  result.park.rmse[,3] <- unlist(py$error)
  result.park.meancrps[,3] <- unlist(py$crps)
  result.park.comptime[,3] <- unlist(py$ctime)
    
  source("Borehole simulation.R")
  py <- py_run_file("python code/Borehole.py")
  result.borehole.rmse[,3] <- unlist(py$error)
  result.borehole.meancrps[,3] <- unlist(py$crps)
  result.borehole.comptime[,3] <- unlist(py$ctime)
    
  source("Currin simulation.R")
  py <- py_run_file("python code/Currin.py")
  result.currin.rmse[,3] <- unlist(py$error)
  result.currin.meancrps[,3] <- unlist(py$crps)
  result.currin.comptime[,3] <- unlist(py$ctime)
    
  source("Franke simulation.R")
  py <- py_run_file("python code/Franke.py")
  result.franke.rmse[,3] <- unlist(py$error)
  result.franke.meancrps[,3] <- unlist(py$crps)
  result.franke.comptime[,3] <- unlist(py$ctime)
  result.franke.rmse[,3][is.nan(result.franke.rmse[,3])] <- NA
  result.franke.meancrps[,3][is.nan(result.franke.meancrps[,3])] <- NA
  result.franke.comptime[,3][is.nan(result.franke.comptime[,3])] <- NA
}else{
  source("Perdikaris simulation.R")
  source("Branin simulation.R")
  source("Park simulation.R")
  source("Borehole simulation.R")
  source("Currin simulation.R")
  source("Franke simulation.R")
}

source("plot Fig7.R")
numerical.rmse
```

<img src="figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

##### Section 5.1: Reproducing Figure S14

This reproduces the CRPS results for synthetic functions.

``` r
source("plot FigS14.R")
numerical.crps
```

<img src="figure-gfm/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

##### Section 5.1: Reproducing Figure 8

This reproduces the computation time results (sec.) for synthetic
functions.

``` r
source("plot Fig8.R")
numerical.time
```

<img src="figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

## Section 5.2

##### Section 5.2: Reproducing Figure 9

This section reproduces the active learning results for the Perdikaris
function. Running ALC can take a significant amount of time. We
recommend setting **ALC=FALSE** during your initial attempt to avoid
long running times.

``` r
# Run 6 active learning strategies for Perdikaris function
source("ALD perdikaris.R")
source("ALM perdikaris.R")
if(ALC) source("ALC perdikaris.R")
source("ALMC perdikaris.R")
source("AL Cokm perdikaris.R")
source("AL MRSUR perdikaris.R")

source("plot Fig9.R")
plotalperd
```

<img src="figure-gfm/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure 10

This reproduces the final RMSE and proportion of low-fidelity data for
Perdikaris function.

``` r
source("plot Fig10.R")
plotalperd2
```

<img src="figure-gfm/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure S15

This section reproduces the active learning results for the Park
function. Running ALC can take a significant amount of time. We
recommend setting **ALC=FALSE** during your initial attempt to avoid
long running times.

``` r
# Run 6 active learning strategies for Park function
source("ALD park.R")
source("ALM park.R")
if(ALC) source("ALC park.R")
source("ALMC park.R")
source("AL Cokm park.R")
source("AL MRSUR park.R")

source("plot FigS15.R")
plotalpark
```

<img src="figure-gfm/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

##### Section 5.2: Reproducing Figure S16

This reproduces the final RMSE and proportion of low-fidelity data for
Park function.

``` r
source("plot FigS16.R")
plotalpark2
```

<img src="figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

## Section 6

##### Section 6: Reproducing Figure 11

This produces the results for for thermal stress analysis of jet engine
turbine blade application. This section includes MATLAB code for running
the finite element simulations for simulating the turbine blade. Before
running .m files, you may need to set a proper working environment.

``` r
# Run Blade simulations
if(NARGP){
  source("Blade simulation.R")
  py <- py_run_file("python code/Blade.py")
  result.blade.rmse[,3] <- unlist(py$error)
  result.blade.meancrps[,3] <- unlist(py$crps)
  result.blade.comptime[,3] <- unlist(py$ctime)
}else{
  source("Blade simulation.R")
}

source("plot Fig11.R")
plot.blade
```

<img src="figure-gfm/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

##### Section 6: Reproducing Figure 12

This section reproduces the active learning results for the Blade data.
Running ALC can take a significant amount of time. We recommend setting
**ALC=FALSE** during your initial attempt to avoid long running times.

``` r
# Run 6 active learning strategies for Blade data
source("ALD Blade.R")
source("ALM Blade.R")
if(ALC) source("ALC Blade.R")
source("ALMC Blade.R")
source("AL Cokm Blade.R")
source("AL MRSUR Blade.R")

source("plot Fig12.R")
plotalblade
```

<img src="figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

##### Section 6: Reproducing Figure 13

This reproduces the plots of final RMSE and proportion of low-fidelity
data for Blade data.

``` r
source("plot Fig13.R")
plotalblade2
```

<img src="figure-gfm/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />
