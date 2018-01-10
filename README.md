# CNORprob

This is a probabilistic logic variant of [CellNOpt](https://www.bioconductor.org/packages/release/bioc/html/CNORode.html) which allows for quantitative optimisation of logical network for (quasi-)steady-state data. The core optimisation pipeline is derived from [FALCON](https://github.com/sysbiolux/FALCON): a toolbox for the fast contextualization of logical networks.

CNORprob takes prior knowledge network and experimental data in SIF and MIDAS formats, respectively, and it outputs graphical representation of results in CellNOpt's format. Post-optimisation analyses including systematic edge-knockout, systematic node-knockout, and local parameter sensitivity analyses (LPSA) were also included in the pipeline as offered in the original FALCON toolbox. 

## Getting Started

A vignette of CNORprob is provided in the 'vignette' folder which could either be used as a Driver script to run the pipeline step-by-step or to simply generate the whole results with Knit(r). A condensed set of optmisation settings are listed at the beginning of the vignette for users' convenience.

### Prerequisites

CNORprob requires prior installation of the "CellNOptR" package (install either from Bioconductor or GitHub page via devtools). It also requires the "Rsolnp" packages as the optimisation algorithm as well as "R.utils" to call for some related functions. Please run the following commands at the beginning of the vignette to perform the installation.

```R
# Install required packages (if were not previously installed)
source("https://bioconductor.org/biocLite.R")
biocLite("CellNOptR") # for loading network/data and pre-processing
# --- OR ---
install.packages("devtools")
library(devtools)
install_github('saezlab/CellNOptR/packages/CellNOptR', force=TRUE)

install.packages("Rsolnp") # optimisation algorithm
install.packages("R.utils") # Timeout control

```

### Installing

CNORprob is currently available for the installation as an R-package from our GitHub page

```R
# Install CNORprob from Github (or load library for development version)
library(devtools)
install_github('saezlab/CellNOptR/packages/CNORprob', force=TRUE)
# --- OR --- #
library(devtools)
load_all()
```

## Running the tests

Several examples are available as the test case for CNORprob. Users can select the examples by assigning the example number to the "Selected_Model" variable. Current examples include: 
1) CNOToy (the running example of the CellNOpt package)
2) FALCON pipeline example, see [paper](https://academic.oup.com/bioinformatics/article/33/21/3431/3897376)
3) PDGF (dissecting PDGF signalling in Gastrointestinal Stromal Tumour - GIST), see [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156223)
4) L-plastin (investigating L-plastin activating signalling pathways in breast cancer cell lines), see [paper](http://www.fasebj.org/content/30/3/1218.long)
5) (Under investigation) Stress-Response (SR) crosstalk between oxidative stress and DNA damage response pathways
6) (Under investigation) CCl4-induced liver injury in mice (in vivo)
7) (Under investigation) APAP-induced Jnk-Jun crosstalk in hepatocyte (in vitro)
8) An example of how to model network motif with XOR gate 

```R
# Select model
# 1 = CNOToy, 2 = FALCON Ex, 3 = PDGF, 4 = L-plastin
# 5 = SR-crosstalk, 6 = CCl4 IFADO, 7 = APAP-JNK, 8 = pp-Raf XOR
Selected_Model <- 1
```

### Setting and options for optimisation and subsequent analyses 

A set of optimisation parameters can be assigned at the beginning of the vignette. Several options can be chosen for the printing of intermedate results and the extraction of the final results. In the post optimisation analyses section, users can define which types of analysis, how many rounds (as well as the number of parameters' increment for LPSA) to perform.

```R
# Assign optimisation and parameter settings
optRound_optim    <- 1        # rounds of optimisation
L1Reg             <- 0.01     # assign weight for L1-regularisation
MaxTime           <- 180      # time for each round of optimisation [seconds]
HLbound           <- 0.5      # cut-off for high and low weights
SSthresh          <- 2e-16    # cut-off for states' difference at steady-state
printCost         <- 0        # print or not print intermediate fitting cost [0,1]
PlotIterations    <- 1        # rounds of optimisation to generate plots
SaveOptResults    <- TRUE     # [TRUE,FALSE] generate reports from optimisation

# Optimiser (rsolnp) options/control list (see rsolnp vignette: https://cran.r-project.org/web/packages/Rsolnp/Rsolnp.pdf)
rho               <- 1        # penalty weighting scaler / default = 1
outer.iter        <- 30       # maximum major iterations / default = 400
inner.iter        <- 30       # maximum minor iterations / default = 800
delta             <- 1e-7     # relative step size eval. / default = 1e-7
tol               <- 1e-8     # relative tol. for optim. / default = 1e-8
trace             <- 1        # print objfunc every iter / default = 1

# Post-optimisation analysis
Analyses          <- c(T,T,T) # [F,T] edge knockout, node knockout, sensitivity analysis
optRound_analysis <- 1        # rounds of optmisation in each analysis
LPSA_Increments   <- 2        # number of increments in LPSA analysis
```

Subsequent to these assignments, users can further run the rest of the script (vignette) on R-console or simply use Knit to generate a HTML report from the pipeline. Note that the results from the optimisation and subsequent analyses will be saved in the sub-folder "Results" within the "vignette" directory

### Additonal/Special assignments in CNORprob

Apart from the core setting options above, users can also assign additional features e.g. the preprocessing of prior knowledge network (please refer to the documentation of CellNOpt for more detail) as well as the network constraint.

"CNORprob_buildModel" provide 4 additional variables: 
1) "expandOR" refers to the expansion of the OR gate which it was not assigned from the initial list. 
2) "ORlist" refers to the introduction of specific OR relationship for specific interaction e.g. 'Grb2SoS_OR_GabSOS=GGSOS' for PDGF model. 
3) "HardConstraint" refers to the assignment that the sum of all weights/probabilities for all incoming activation reactions to be 1 (might be too strict for some cases). 
4) "Force" refers to the assignment of weight/probability for a single activated to always be 1 (might also be too strict for several cases).

```R
# Preprocessing (please see CellNopt documentation)
model <- preprocessing(CNOlist, pknmodel, expansion=FALSE,compression=TRUE, cutNONC=TRUE, verbose=FALSE)

# Model building with specified OR gate interaction
estim <- CNORprob_buildModel(CNOlist,model,expandOR=FALSE,ORlist=c("Grb2SOS_OR_GabSOS=GGSOS"),HardConstraint=TRUE,Force=TRUE,L1Reg=L1Reg,HLbound=HLbound,SSthresh=SSthresh,PlotIterations=PlotIterations,rsolnp_options=rsolnp_options)
```

## Authors

Panuwat Trairatphisan (panuwat.trairatphisan -at - gmail.com)

See also the list of [contributors](https://github.com/saezlab/CNORprob/contributors) who participated in this project.

## Acknowledgement

CNORprob have been developed for the modelling projects within the [TransQST Consortium](https://transqst.org)

"This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant agreement No 116030. The Joint Undertaking receives support from the European Union's Horizon 2020 research and innovation programme and EFPIA."
