# CNORprob

This is a probabilistic logic variant of [CellNOpt](https://www.bioconductor.org/packages/release/bioc/html/CNORode.html) which allows for quantitative optimisation of logical network for (quasi-)steady-state data. The core optimisation pipeline is derived from [FALCON](https://github.com/sysbiolux/FALCON): a toolbox for the fast contextualization of logical networks.

CNORprob takes prior knowledge network and experimental data in SIF and MIDAS formats, respectively, and it outputs graphical representation of results in CellNOpt's format. Post-optimisation analyses including systematic edge-knockout, systematic node-knockout, bootstrapping and local parameter sensitivity analyses (LPSA) were also included in the pipeline as offered in the original FALCON toolbox. 

## Getting Started

A vignette of CNORprob is provided in the 'vignette' folder which could either be used as a Driver script to run the pipeline step-by-step or to simply generate the whole results with Knit(r). A condensed set of optmisation settings are listed at the beginning of the vignette for users' convenience.

### Prerequisites

CNORprob requires prior installation of the "CellNOptR" package (install either from Bioconductor or GitHub page via devtools). It also requires the "Rsolnp" packages as the optimisation algorithm as well as "R.utils" to call for some related functions. Please run the following commands at the beginning of the vignette to perform the installation.

```R
# Install required packages (if were not previously installed)
# CellNOptR
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CellNOptR", version = "3.8")
# --- OR ---
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/CellNOptR', force=TRUE)

# Rsolnp and R.utils
install.packages("Rsolnp") # optimisation algorithm
install.packages("R.utils") # Timeout control
```

### Installing

CNORprob is currently available for the installation as an R-package from our GitHub page

```R
# Install CNORprob from Github (or load library for development version)
# install.packages('devtools') # in case devtools hasn't been installed
library(devtools)
install_github('saezlab/CNORprob') # Doesn't work currently as the repository is still private
# or download the source file from GitHub and install from source
install.packages('path_to_extracted_CNORprob_directory', repos = NULL, type="source")
```

## Running the tests

Several examples are available as the test case for CNORprob. Users can select the examples by assigning the example number to the "Selected_Model" variable. Current examples include: 
1) CNOToy (the running example of the CellNOpt package)
2) FALCON pipeline example, see [paper](https://academic.oup.com/bioinformatics/article/33/21/3431/3897376)
3) PDGF (dissecting PDGF signalling in Gastrointestinal Stromal Tumour - GIST), see [paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156223)
4) L-plastin (investigating L-plastin activating signalling pathways in breast cancer cell lines), see [paper](http://www.fasebj.org/content/30/3/1218.long)

```R
# Load required packages
library(CNORprob)
library(CellNOptR) # for loading network/data and pre-processing
library(Rsolnp) # optimisation algorithm
library(R.utils) # timeout control

# Select CNORprob example
# 1 = CNO-Toy, 2 = FALCON-Ex, 3 = PDGF, 4 = L-plastin
CNORprob_Example <- 1
CNORprob_inputs <- loadExample_Prob(CNORprob_Example = CNORprob_Example)
optmodel <- CNORprob_inputs$optmodel
optCNOlist <- CNORprob_inputs$optCNOlist
```

### Setting and options for optimisation and subsequent analyses 

A set of optimisation parameters can be assigned at the beginning of the vignette. Several options can be chosen for the printing of intermedate results and the extraction of the final results. In the post optimisation analyses section, users can define which types of analysis, how many rounds (as well as the number of parameters' increment for LPSA) to perform.

```R
# Assign optimisation and parameter settings
optRound_optim    <- 3        # rounds of optimisation
L1Reg             <- 1e-4     # assign weight for L1-regularisation
MaxTime           <- 180      # time for each round of optimisation [seconds]
HLbound           <- 0.5      # cut-off for high and low weights
SSthresh          <- 2e-16    # cut-off for states' difference at steady-state
printCost         <- 0        # print or not print intermediate fitting cost [0,1]
PlotIterations    <- 1        # rounds of optimisation to generate plots
SaveOptResults    <- TRUE     # [TRUE,FALSE] generate reports from optimisation

# Optimiser (rsolnp) options/control list (see rsolnp vignette: https://cran.r-project.org/web/packages/Rsolnp/Rsolnp.pdf)
rho               <- 1        # penalty weighting scaler / default = 1
outer.iter        <- 100       # maximum major iterations / default = 400
inner.iter        <- 100       # maximum minor iterations / default = 800
delta             <- 1e-7     # relative step size eval. / default = 1e-7
tol               <- 1e-8     # relative tol. for optim. / default = 1e-8
trace             <- 1        # print objfunc every iter / default = 1

# Post-optimisation analysis
Analyses          <- c(T,T,T,T) # [F,T] edge knockout, node knockout, sensitivity analysis
optRound_analysis <- 1        # rounds of optmisation in each analysis
LPSA_Increments   <- 2        # number of increments in LPSA analysis
BS_Type           <- 1        # Type of Bootstrapping [1=resample with replacement from residual; 2=resampling from mean & variant]
BS_Round          <- 5       # number of rounds for bootstrapping
```

Subsequent to these assignments, users can run the vignette using Knit to generate a HTML report from the pipeline. In this case, the results from the optimisation and subsequent analyses will be saved in the sub-folder "Results" within the "vignette" directory. Alternatively, users can also proceed further with manual execution of the pipeline on R-console as outlined below.

### Optimisation

Prior to the optimisation, model and experimental descriptions together with optimisation parameters are combined into an optimisation object 'estim', then the optimisation can be performed.

```R
estim_Result   <<- list() # Initialise global variable of results
rsolnp_options  <- list(rho=rho,outer.iter=outer.iter,inner.iter=inner.iter,delta=delta,tol=tol,trace=trace) # Collapase optimisation options

# Generate optimisation object
estim <- CNORprob_buildModel(optCNOlist,optmodel,expandOR=CNORprob_inputs$ProbExpandOR,ORlist=CNORprob_inputs$ORlist,HardConstraint=CNORprob_inputs$ProbHardConstraint,Force=CNORprob_inputs$ProbForce,L1Reg=L1Reg,HLbound=HLbound,SSthresh=SSthresh,PlotIterations=PlotIterations,rsolnp_options=rsolnp_options)

# Run Optimisation
estim$maxtime <- MaxTime; estim$printCost <- printCost
estim$optimOptions <- c(rho,outer.iter,inner.iter,delta,tol,trace)
res <- CNORprob_optimise(estim,optRound_optim,SaveOptResults)
```

### Plot results and post-optimisation analyses

The results from CNORprob are stored in the 'res' variable and can be used to plot resulting figures (compare fitting quality and fitted network with weights on edges).

```R
# Plot results
CNORprob_plotFit(optmodel,optCNOlist,estim,res,show=TRUE, plotPDF=TRUE, tag=NULL, plotParams=list(cex=0.8, cmap_scale=1, ymin=0))
MappedProb <- CNORprob_mapModel(optmodel,optCNOlist,estim,res)
plotModel(MappedProb$model,MappedProb$CNOlist,round(MappedProb$bString,digits=2))
```

In addition, the 'res' variable can also be passed to run post-optimisation analyses i.e. local parameter sensitivity analysis (LPSA), edge/edge knockout and boot-strapping analysis. Results of plotting and post-optimisation analyses will be stored in the "Results" folder.

```R
# Post-optimisation analyses
estim$ProbCompression <- CNORprob_inputs$ProbCompression
estim$ProbCutNONC <- CNORprob_inputs$ProbCutNONC
estim$ProbExpandOR <- CNORprob_inputs$ProbExpandOR
estim$optRound_analysis <- optRound_analysis
estim_original <- estim
estim_based <- estim_original; if (Analyses[1]) { estim_Result  <- CNORprob_edgeKO(optmodel,optCNOlist,estim_based,res) }
estim_based <- estim_original; if (Analyses[2]) { estim_Result  <- CNORprob_nodeKO(optmodel,optCNOlist,estim_based,res) }
estim_based <- estim_original; if (Analyses[3]) { estim_Result  <- CNORprob_LPSA(estim_based,res,HLbound,LPSA_Increments,Force=F) }
estim_based <- estim_original; if (Analyses[4]) { estim_Result  <- CNORprob_BS(optmodel,optCNOlist,estim_based,res,BS_Type,BS_Round) }

save(estim_Result,file="Results/CNORprob_PostHocResults.Rdata")
```

### Manual assignment of input files & additonal/special assignments in CNORprob

In case users would like to use own prior knowledge network (PKN) and experimental data which are in the MIDAS format (please refer to the original documentation of the CellNOptR package and publication), users could the following codes to manually assign the input files which will then converted into model and data in the CellNOpt format.

```R
# Manual assignment of network (PKN) and experimental data (MIDAS) files
pknmodel <- readSIF('path_to_your_PKN_model_here') # build a prior-knowledge network from SIF file
CNOlist <- CNOlist('path_to_your_MIDAS_data_here') # import experimental data from MIDAS file
``` 

Apart from applyaing the default setting as used in CNORprob examples, users can also assign additional features e.g. the preprocessing of prior knowledge network (please refer to the documentation of CellNOpt for more detail) as well as the network constraint.

"CNORprob_buildModel" provide 4 additional variables: 
1) "expandOR" refers to the expansion of the OR gate which it was not assigned from the initial list. 
2) "ORlist" refers to the introduction of specific OR relationship for specific interaction e.g. 'Grb2SoS_OR_GabSOS=GGSOS' for PDGF model. 
3) "HardConstraint" refers to the assignment that the sum of all weights/probabilities for all incoming activation reactions to be 1 (might be too strict for some cases). 
4) "Force" refers to the assignment of weight/probability for a single activated to always be 1 (might also be too strict for several cases).

```R
# Assign CNORprob-specitifc parameters for pre-processing step
ProbCompression <- F;ProbCutNONC <- F; ProbExpandOR <- F; ProbHardConstraint <- F; ProbForce <- F # Default (most-relaxed) setting 

# Run CNORprob-specific pre-processing (expansion=FALSE by default and report)
ModDatppProb <- preprocessing_Prob(CNOlist, pknmodel, expansion=FALSE,
                                     compression=ProbCompression, cutNONC=ProbCutNONC,
                                     verbose=FALSE) 
optmodel <- ModDatppProb$cutModel
optCNOlist <- ModDatppProb$data
``` 

The variables 'optmodel' and 'optCNOlist' could then be passed into the optimisation pipeline to generate results as outlined above (start from the 'Optimisation' section).

## Authors

Panuwat Trairatphisan (panuwat.trairatphisan -at - gmail.com)

See also the list of [contributors](https://github.com/saezlab/CNORprob/contributors) who participated in this project.

## License

Distributed under the GNU GPLv3 License. See accompanying file [LICENSE.txt](https://github.com/saezlab/combiMS/blob/master/LICENSE.txt) or copy at [http://www.gnu.org/licenses/gpl-3.0.html](http://www.gnu.org/licenses/gpl-3.0.html).

## References

[De Landtsheer et al.](https://academic.oup.com/bioinformatics/article/33/21/3431/3897376):

> De Landtsheer S, Trairatphisan P, Lucarelli P, and Sauter, T. (2017). FALCON: a toolbox for the fast contextualization of logical networksa. *Bioinformatics*, Volume 33, Issue 21, 1 November 2017, Pages 3431–3436, https://doi.org/10.1093/bioinformatics/btx380.

## Acknowledgement

CNORprob have been developed for the modelling projects within the [TransQST Consortium](https://transqst.org)

"This project has received funding from the Innovative Medicines Initiative 2 Joint Undertaking under grant agreement No 116030. The Joint Undertaking receives support from the European Union's Horizon 2020 research and innovation programme and EFPIA."
