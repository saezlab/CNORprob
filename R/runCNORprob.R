#' @title runCNORprob
#'
#' @description A one-linear wrapper of the CNORprob pipeline
#'
#' @param model A filename of prior knowledge network (PKN) in the SIF format
#' @param data A measurement filename in the MIDAS format
#' @param CNORprob_Example A variable to select CNORprob examples (1-4)
#' @param ProbCompression A variable to define whether to perform a compression step 
#' @param ProbCutNONC A variable to define whether to cut non-observable non-controllable node from PKN
#' @param ProbExpandOR A variable to define whether to expand OR gate in the CNOR format
#' @param ORlist A list of interaction with OR gate (leave NULL is none)
#' @param ProbHardConstraint Apply hard constraints on probability parameters (sum of all activatory edges = 1)
#' @param ProbForce A variable to define whether to force all single activatory edge to a node to always be 1
#' @param optRound_optim Number of rounds for optimisation
#' @param L1Reg Weight for L1-regularisation
#' @param MaxTime Time for each round of optimisation (seconds)
#' @param HLbound A cut-off for high and low weights
#' @param SSthresh A cut-off for states' difference to define whether steady-state is reached
#' @param printCost A variable to select whether to print or not to print intermediate fitting cost (0,1)
#' @param PlotIterations Number of rounds for optimisation to generate plots
#' @param SaveOptResults A variable to define whether to generate reports from optimisation (TRUE,FALSE)
#' @param rho Penalty weighting scaler / default = 1
#' @param outer.iter Maximum major iterations / default = 400
#' @param inner.iter Maximum minor iterations / default = 800
#' @param delta Relative step size eval. / default = 1e-7
#' @param tol Relative tol. for optim. / default = 1e-8
#' @param iter Print objfunc every itereration / default = 1
#' @param Analysis_edgeKO Perform edge KO analysis (T/F)
#' @param Analysis_nodeKO Perform node KO analysis (T/F)
#' @param Analysis_LPSA Perform local parameter sensitivity analysis (T/F)
#' @param Analysis_BS Perform bootstrapping analysis (T/F)
#' @param optRound_analysis Number of rounds for optmisation in each post-optimisation analysis
#' @param LPSA_Increments Number of increments in LPSA analysis
#' @param BS_Type Type of Bootstrapping (1=resample with replacement from residual; 2=resampling from mean & variant)
#' @param BS_Round Number of rounds for bootstrapping
#' 
#' 
#' @export

runCNORprob = function(model=NULL,
                       data=NULL,
                       CNORprob_Example=1,
                       ProbCompression=F,
                       ProbCutNONC=F,
                       ProbExpandOR=F,     
                       ORlist=NULL,
                       ProbHardConstraint=F,      
                       ProbForce=F,
                       optRound_optim=3,        
                       L1Reg=1e-4,     
                       MaxTime=180,      
                       HLbound=0.5,      
                       SSthresh=2e-16,   
                       printCost=0,        
                       PlotIterations=1,        
                       SaveOptResults=TRUE,     
                       rho=1,        
                       outer.iter=100,       
                       inner.iter=100,       
                       delta=1e-7,     
                       tol=1e-8,     
                       trace=1,        
                       Analysis_edgeKO=T,
                       Analysis_nodeKO=T,
                       Analysis_LPSA=T,
                       Analysis_BS=T,
                       optRound_analysis=1,       
                       LPSA_Increments=2,        
                       BS_Type=1,        
                       BS_Round=5       
                       ) {

  library(CNORprob)
  library(CellNOptR) # for loading network/data and pre-processing
  library(Rsolnp) # optimisation algorithm
  library(R.utils) # timeout control
  
  # Select CNORprob example
  # 1 = CNO-Toy, 2 = FALCON-Ex, 3 = PDGF, 4 = L-plastin
  # CNORprob_Example <- 1
  if (!is.null(CNORprob_Example)) {
    CNORprob_inputs <- loadExample_Prob(CNORprob_Example = CNORprob_Example)
    optmodel <- CNORprob_inputs$optmodel
    optCNOlist <- CNORprob_inputs$optCNOlist
  } else {
    if (!is.null(model) & !is.null(data)) {
      # Manual assignment of network (PKN) and experimental data (MIDAS) files
      pknmodel <- readSIF(model) # build a prior-knowledge network from SIF file
      CNOlist <- CNOlist(data) # import experimental data from MIDAS file
      # Now use CNORprob specific pre-processing (expansion=FALSE by default and report)
      ModDatppProb <- preprocessing_Prob(CNOlist, pknmodel, expansion=FALSE, # model in CNORprob shouldn't be expanded by pre-processing step
                                         compression=ProbCompression, cutNONC=ProbCutNONC,
                                         verbose=FALSE) 
      optmodel <- ModDatppProb$cutModel
      optCNOlist <- ModDatppProb$data
      CNORprob_inputs <- NULL
      CNORprob_inputs$optmodel=optmodel
      CNORprob_inputs$optCNOlist=optCNOlist
      CNORprob_inputs$ProbExpandOR=ProbExpandOR
      CNORprob_inputs$ProbHardConstraint=ProbHardConstraint
      CNORprob_inputs$ProbForce=ProbForce
      CNORprob_inputs$ORlist=ORlist
    } else {
      stop("Please provide the prior knowledge network (PKN) file and experimental data in MIDAS format or select a CNORprob example")
    }
  }

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
  pdf("Results/Figures/Optimised_CNORprob_Model.pdf")
  plotModel(MappedProb$model,MappedProb$CNOlist,round(MappedProb$bString,digits=2))
  dev.off()
  
  # Post-optimisation analyses
  estim$ProbCompression <- CNORprob_inputs$ProbCompression
  estim$ProbCutNONC <- CNORprob_inputs$ProbCutNONC
  estim$ProbExpandOR <- CNORprob_inputs$ProbExpandOR
  estim$optRound_analysis <- optRound_analysis
  estim_original <- estim
  estim_based <- estim_original; if (Analysis_edgeKO) { estim_Result  <- CNORprob_edgeKO(optmodel,optCNOlist,estim_based,res) }
  estim_based <- estim_original; if (Analysis_nodeKO) { estim_Result  <- CNORprob_nodeKO(optmodel,optCNOlist,estim_based,res) }
  estim_based <- estim_original; if (Analysis_LPSA) { estim_Result  <- CNORprob_LPSA(estim_based,res,HLbound,LPSA_Increments,Force=F) }
  estim_based <- estim_original; if (Analysis_BS) { estim_Result  <- CNORprob_BS(optmodel,optCNOlist,estim_based,res,BS_Type,BS_Round) }
  
  save(estim_Result,file="Results/CNORprob_PostHocResults.Rdata")
  
  return(estim_Result)

}

# ======= End of the script ====== #
