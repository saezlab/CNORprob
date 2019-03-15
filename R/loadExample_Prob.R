#' @title loadExample_Prob
#'
#' @description loadExample_Prob
#'
#' @export

#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EMBL - European Bioinformatics Institute
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################

loadExample_Prob<-function(CNORprob_Example=NULL){

  if (!is.numeric(CNORprob_Example)){stop("Please provide a number of CARNIVAL example [1-4]")}

  # Load required libraries
  library(CellNOptR)
  library(CNORprob)
  
  Selected_Model <- CNORprob_Example 
  
  CNORprob_Example_Names <- c("CNO-Toy","FALCON-Ex","PDGF","L-plastin")
  
  # Assign selected model to pknmodel and CNOlist
  if (Selected_Model==1) {
    model <- system.file("ToyPKNMMB.sif",package="CNORprob")
    data  <- system.file("ToyDataMMB.csv",package="CNORprob")
  } else if (Selected_Model==2) {
    model <- system.file("FALCON_Ex1_Model.sif",package="CNORprob")
    data  <- system.file("FALCON_Ex1_Exp.csv",package="CNORprob")
  } else if (Selected_Model==3) {
    model <- system.file("FALCON_PDGF_Model.sif",package="CNORprob")
    data  <- system.file("FALCON_PDGF_Exp.csv",package="CNORprob")
  } else if (Selected_Model==4) {
    model <- system.file("FALCON_LPL_BT20_Model.sif",package="CNORprob")
    data  <- system.file("FALCON_LPL_BT20_Exp.csv",package="CNORprob")
  } else {
    print('::: Please select the model number from the list :::')
  }
  
  pknmodel <- readSIF(model) # build a prior-knowledge network from SIF file
  CNOlist <- CNOlist(data) # build a CNOlist from data file in MIDAS format
  
  # pre-processing based on the selected model
  if (Selected_Model==1) {ProbCompression <- T;ProbCutNONC <- T; ProbExpandOR <- T; ProbHardConstraint <- F; ProbForce <- F
  } else {ProbCompression <- F;ProbCutNONC <- F; ProbExpandOR <- F; ProbHardConstraint <- F; ProbForce <- F}
  
  if (Selected_Model==3) {ORlist=c("Grb2SOS_OR_GabSOS=GGSOS"); # Manually assign OR gate
    tol <- 1e-4 # reduced tolerance of solution to speed up optimisation time
  } else if (Selected_Model==4) {tol <- 1e-4; ORlist=NULL} else {tol <- 1e-8; ORlist=NULL}
  
  # Now use CNORprob specific pre-processing (expansion=FALSE by default and report)
  ModDatppProb <- preprocessing_Prob(CNOlist, pknmodel, expansion=FALSE,
                                     compression=ProbCompression, cutNONC=ProbCutNONC,
                                     verbose=FALSE) 
  optmodel <- ModDatppProb$cutModel
  optCNOlist <- ModDatppProb$data
  
  print(paste0("CNORprob example Number ",CNORprob_Example," - ",CNORprob_Example_Names[CNORprob_Example]," - is loaded."))
  print(paste0("Please check the model and experimental descriptions in the variables 'optmodel' and 'optCNOlist' respectively"))
              
  return(list(optmodel=optmodel,optCNOlist=optCNOlist,ProbExpandOR=ProbExpandOR,ProbHardConstraint=ProbHardConstraint,ProbForce=ProbForce,ORlist=ORlist,tol=tol))
}
