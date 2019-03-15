## ---- message=FALSE, warning=FALSE---------------------------------------

rm(list=ls()); cat("\014") # clean screen and variables 

# Load prior library and R-scripts
library(CNORprob)
library(CellNOptR) # for loading network/data and pre-processing
library(Rsolnp) # optimisation algorithm
library(R.utils) # Timeout control

# Select CNORprob example
# 1 = CNO-Toy, 2 = FALCON-Ex, 3 = PDGF, 4 = L-plastin
CNORprob_Example <- 1
CNORprob_inputs <- loadExample_Prob(CNORprob_Example = CNORprob_Example)
optmodel <- CNORprob_inputs$optmodel
optCNOlist <- CNORprob_inputs$optCNOlist

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

# ================================================ # 
# ======= Run the driver with option+cmd+R ======= #
# ================================================ # 

estim_Result   <<- list() # Initialise global variable of results
rsolnp_options  <- list(rho=rho,outer.iter=outer.iter,inner.iter=inner.iter,delta=delta,tol=tol,trace=trace) # Collapase optimisation options
  
# Generate optimisation object
estim <- CNORprob_buildModel(optCNOlist,optmodel,expandOR=CNORprob_inputs$ProbExpandOR,ORlist=CNORprob_inputs$ORlist,HardConstraint=CNORprob_inputs$ProbHardConstraint,Force=CNORprob_inputs$ProbForce,L1Reg=L1Reg,HLbound=HLbound,SSthresh=SSthresh,PlotIterations=PlotIterations,rsolnp_options=rsolnp_options)

# Run Optimisation
estim$maxtime <- MaxTime; estim$printCost <- printCost
estim$optimOptions <- c(rho,outer.iter,inner.iter,delta,tol,trace)
res <- CNORprob_optimise(estim,optRound_optim,SaveOptResults)

# Plot results
CNORprob_plotFit(optmodel,optCNOlist,estim,res,show=TRUE, plotPDF=TRUE, tag=NULL, plotParams=list(cex=0.8, cmap_scale=1, ymin=0))
MappedProb <- CNORprob_mapModel(optmodel,optCNOlist,estim,res)
plotModel(MappedProb$model,MappedProb$CNOlist,round(MappedProb$bString,digits=2))

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

# ===== End of the script ===== #


