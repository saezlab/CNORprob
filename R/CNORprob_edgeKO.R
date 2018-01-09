#' @title CNORprob_edgeKO
#'
#' @description Systematically knock-out each single edge from the network, run optimisation of each network variant and compare models with fitting costs and AIC
#'
#' @export

CNORprob_edgeKO = function(model,CNOlist,estim,res) {

  # estim_orig <<- estim
  estim_orig <- estim
  model_orig <- model # Keep original CNOR model

  optRound_KO <- estim$optRound_analysis
  bestx <- res$BestFitParams
  num_params <- length(estim$param_vector)
  bestcost <- res$OptResults[1,1]
  L1Reg <- estim$L1Reg
  SSthresh <- estim$SSthresh

  # AIC calculation
  N <- length(estim$Output_vector)
  SSE <- bestcost - (L1Reg*sum(res$ParamResults[,2]!=0))
  SSE_orig <- SSE
  p <- num_params

  AIC_complete <- N*log(SSE/N) + 2*p # AIC for base model

  NrReacID <- length(model$reacID)

  p_KD <- rep(0,1,NrReacID)

  cost_KD <- rep(0,1,NrReacID)
  AIC_KD <- rep(0,1,NrReacID)

  N_r_All <- NULL
  p_r_All <- NULL
  res <- NULL

  # pb <- tkProgressBar(title = "Edge Knock-out analysis", min = 0, max = length(p_KD), width=300)

  for (counter in  1:length(p_KD)) {

    print("=================================")
    print(paste("Knocking down parameter:", toString(counter), "/", toString(length(p_KD))))
    print("=================================")

    print(model_orig$reacID[counter])

    model_KD <- model_orig

    # setTkProgressBar(pb, counter, label=paste( round(counter/length(p_KD)*100, 0), "% done"))
    # setTkProgressBar(pb, counter, label=paste("Knocking down parameter:", toString(counter), "/", toString(length(p_KD))))

    model_KD$reacID <- model$reacID[-counter]
    Target_KOed <- strsplit(x = model$reacID[counter],split = "=")[[1]][2]

    if (model$reacID[1] == "EGF=PI3K") { # if CNOToy model
      estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=TRUE,HardConstraint=FALSE,Force=FALSE)
    } else if (model$reacID[1] == "PDGFL=PDGFR") { # if PDGF modeel
      # estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE,ORlist=c("Grb2SOS_OR_GabSOS=GGSOS"))
      estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE)
    } else {
      estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE)
    }

    # If the KOed interaction left with target being not activated anymore -> Fix target to zero (instead of random number)
    if (model$reacID[1] == "EGF=PI3K") { # if CNOToy model with expandOR
      if (sum(grepl(paste("\\b",Target_KOed,"\\b",sep=""),estim$Interactions[,3],fixed = TRUE))==0) {
        Target_KOed_Idx <- which(Target_KOed==estim$state_names)
        estim$Input_vector <- cbind(estim$Input_vector,rep(0,dim(estim$Input_vector)[1]))
        estim$Input_index  <- cbind(estim$Input_index,rep(Target_KOed_Idx,dim(estim$Input_index)[1]))
      }
    } else {
      if (sum(grepl(Target_KOed,estim$Interactions[,3],fixed = TRUE))==0) {
        Target_KOed_Idx <- which(Target_KOed==estim$state_names)
        estim$Input_vector <- cbind(estim$Input_vector,rep(0,dim(estim$Input_vector)[1]))
        estim$Input_index  <- cbind(estim$Input_index,rep(Target_KOed_Idx,dim(estim$Input_index)[1]))
      }
    }
    estim$maxtime <- estim_orig$maxtime
    estim$printCost <- estim_orig$printCost

    res <- NULL

    tryCatch(
    res <- CNORprob_optimise(estim,optRound=optRound_KO)
    ,
    error=function(e) print("Solver failed... continue the next optimisation round"))

    if (!is.null(res)) {

      CNORprob_plotFit(model_KD,CNOlist,estim,res,plotPDF = FALSE)

      SSE_L1 <- res$OptResults[1,1] # fetch the best cost from optimizations
      SSE <- SSE_L1 - (L1Reg*sum(res$ParamResults[,2]!=0))
      cost_KD[counter] <- SSE # fetch the best cost from optimizations

      # reduced model (-1 parameter)

      N_r = length(estim$Output_vector) # number of datapoints
      p_r= length(estim$param_vector) # number of parameters

      N_r_All <- c(N_r_All,N_r)
      p_r_All <- c(p_r_All,p_r)

      AIC_KD[counter] = N_r*log(cost_KD[counter]/N_r) + 2*p_r

    } else { # if optimiser failed
      cost_KD[counter] <- NA
      AIC_KD[counter] <- NA
    }

  }

  # close(pb)

  # Plot AIC values
  estim <- estim_orig

  # par(mfrow=c(1,1))
  All_AIC <- c(AIC_complete,AIC_KD)
  # All_KO_ParamNames <- c("Orig",estim$param_vector)
  All_KO_ParamNames <- c("Orig",model$reacID)
  cols <- c("gray","green")[(All_AIC < All_AIC[1])+1]
  cols[1] <- "blue" # First column always in blue
  barplot(All_AIC, main="Edge-KO analysis", names.arg = All_KO_ParamNames, col = cols, ylab = "AIC",cex.names = 0.7,las=2)
  abline(h=All_AIC[1],col="red")

  All_Cost <- c(SSE_orig,cost_KD)
  # All_KO_ParamNames <- c("Orig",estim$param_vector)
  All_KO_ParamNames <- c("Orig",model$reacID)
  cols <- c("gray","red")[(All_Cost < All_Cost[1])+1]
  cols[1] <- "blue" # First column always in blue
  barplot(All_Cost, main="Edge-KO analysis", names.arg = All_KO_ParamNames, col = cols, ylab = "SSE",cex.names = 0.8,las=2)
  abline(h=All_Cost[1],col="red")

  # Exporting figure
  pdf("Results/Figures/CNORprob_edgeKO_AIC.pdf")
  barplot(All_AIC, main="Edge-KO analysis", names.arg = All_KO_ParamNames, col = cols, ylab = "AIC",cex.names = 0.7,las=2)
  abline(h=All_AIC[1],col="red")
  dev.off()
  pdf("Results/Figures/CNORprob_edgeKO_Cost.pdf")
  barplot(All_Cost, main="Edge-KO analysis", names.arg = All_KO_ParamNames, col = cols, ylab = "SSE",cex.names = 0.8,las=2)
  abline(h=All_Cost[1],col="red")
  dev.off()

  estim_Result$edgeKO$AIC_KD  <- All_AIC
  estim_Result$edgeKO$Cost_KD <- All_Cost
  estim_Result               <<- estim_Result
  estim                      <<- estim_orig
  estim                      <<- estim

  return(estim_Result)

}

# ======= End of the script ======= #
