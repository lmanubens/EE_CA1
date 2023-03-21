# Analysis and modeling of CA1 pyramidal neurons for the article: Deficits in neuronal architecture but not over-inhibition are main determinants of reduced neuronal network activity in a mouse model of overexpression of Dyrk1A


<!--[![biorXiv shield](https://img.shields.io/badge/arXiv-1709.01233-red.svg?style=flat)](https://www.biorxiv.org/content/10.1101/2023.03.09.531874v1) -->


## Contents

- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Issues](https://github.com/lmanubens/EE_CA1/issues)
- [Citation](#citation)

# Repo Contents

- [T2N modeling code](./model): T2N single neuron multi-comparmental simulation code.
- [Shiny app code](./shiny_app): `Shiny` app code.
- [Test data](./data): Reconstruction set for replicating the results of the manuscript. 

# System Requirements

## Hardware Requirements

The Shiny app requires only a standard computer with enough RAM to support the operations defined by a user. The T2N simulation scripts can also be run in a standard computer. For minimal performance, this will be a computer with about 500 MB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 1.5+ GB  
CPU: single core, 1.6+ GHz

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 8 cores@2.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

The package development version is tested on *Linux*, *MacOS* and *Windows* operating systems. The developmental version of the package has been tested on the following systems:

Linux: Centos 8  
Mac OSX:  Monterey 12.3.1  
Windows:  10 and 11

Both the T2N modeling code and the Shiny app should be compatible with Windows, Mac, and Linux operating systems. Some attention may be needed for path specifications in scripts. 


Before setting up the Shiny app, users should have `R` version 4.1.0 or higher, and several packages set up from CRAN. 

#### Installing R version 4.1.0 on Centos 8

the latest version of R can be installed by adding the EPEL repository to `yum`:

```
sudo yum install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
sudo dnf config-manager --set-enabled powertools

sudo bash -c "$(curl -L https://rstd.io/r-install)"

vim ~/.bashrc
export PATH=$PATH:/opt/R/4.1.0/bin/

source .bashrc 
```

#### Installing R on MacOS

Download and install the corresponding binary .pkg file from:

https://cran.r-project.org/bin/macosx/

#### Installing R on Windows

Download and install the corresponding binary .exe file from:

https://cran.r-project.org/bin/windows/base/


In all cases it should install in about 20 seconds.


# Installation Guide

## T2N modeling installation

Users should install NEURON and T2N following the documentation available at: [https://github.com/MarcelBeining/T2N](https://github.com/MarcelBeining/T2N)

## Shiny app installation

### Package dependencies

Users should install the following packages prior to running the Shiny app, from an `R` terminal:

```
options("repos" = c("CRAN" = "https://cran.rstudio.com",
                    "rforge" = "http://R-Forge.R-project.org"))
                    
install.packages(c('shiny', 'plotly', 'ggpubr', 'afex', 'dplyr', 'reshape2', 'xtable', 'emmeans', 'ggsignif', 
                  'plotrix', 'scales', 'rstatix', 'Matrix'))
```

which will install in about 5 minutes on a machine with the recommended specs.

The Shiny app functions with all packages in the following versions:
```
shiny_1.6.0
plotly_4.9.4.1
ggpubr_0.4.0
afex_0.28-1
dplyr_1.0.7
reshape2_1.4.4
xtable_1.8-4
emmeans_1.6.1
ggsignif_0.6.2
plotrix_3.8-2
scales_1.2.1
rstatix_0.7.0
Matrix_1.3-3
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/lmanubens/EE_CA1/issues). 

### Shiny app Installation

First download the repository code and data locally, by clicking "Code Download ZIP" or typing in a command line:

```
git clone https://github.com/lmanubens/EE_CA1.git
```

From an `R` session in the EE_CA1/shiny_app folder, type:

```
library(shiny)
runApp()
```

# Demo

## T2N modeling scripts

### Connectivity repertoire

To obtain connectivity repertoires, run the script 'connectivity_repertoires.m' from the folder EE_CA1/model, the results are saved in a file called "wiring_metrics.csv".
The script takes about 4 seconds.

### Input-output frequencies
To obtain input-output frequency curves, run the command "export MATLAB_SHELL=/bin/bash", or equivalent, from a terminal. Then run Matlab from terminal and run the script 'freqinout_stimregion_BC.m' from the folder EE_CA1/model. Different sections can be run separately to account for synaptic densities or only dendritic morphological features, and for input to Stratum Radiatum alone or both Stratum Radiatum and Lacunosum. The results are saved in a file called "freqinout_stimregion_BC.csv".

Running each set of simulations takes about 20 minutes using 8 parallel workers.

### Power spectrum with inhibitory feedback
To obtain power spectrums from the set of reconstructed neurons with inhibitory feedback, run the command "export MATLAB_SHELL=/bin/bash", or equivalent, from a terminal. Then run Matlab from terminal and run the script 'net_morph_and_inh.m' or 'net_only_inhibition' from the folder EE_CA1/model. In this case both Stratum Radiatum and Lacunosum are stimulated. The results are saved in a file called "pspectrum_x.csv".

Running each set of simulations takes about 3 minutes using a single core.

## Ready-to-use Shiny web app

For interactive usage of the web app please check:

```
https://linusmg.shinyapps.io/EE_TgDyrk1A/
```
The tabs allow users to choose specific sets of analyses of interest.
The checklists in the top of the web app allow the users to choose specific metrics and data subsets of interest respectively.
Each plot should be generated in a few seconds.
