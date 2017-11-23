#' @title CNORprob_nodeKO
#'
#' @description Systematically knock-out each single node from the network, run optimisation of each network variants and compare models with fitting costs and AIC
#'
#' @export

CNORprob_nodeKO = function(model,CNOlist,estim,res) {

  estim_orig <<- estim
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

  # NrReacID <- length(model$reacID) # EdgeKO
  NrReacID <- length(model$namesSpecies) # NodeKO

  p_KD <- rep(0,1,NrReacID)

  cost_KD <- rep(0,1,NrReacID)
  AIC_KD <- rep(0,1,NrReacID)

  N_r_All <- NULL
  p_r_All <- NULL
  res <- NULL

  # pb <- tkProgressBar(title = "Edge Knock-out analysis", min = 0, max = length(p_KD), width=300)

  for (counter in  1:length(p_KD)) {

    print("=================================")
    print(paste("Knocking down node:", toString(counter), "/", toString(length(p_KD))))
    print("=================================")

    model_KD <- model_orig

    # setTkProgressBar(pb, counter, label=paste( round(counter/length(p_KD)*100, 0), "% done"))
    # setTkProgressBar(pb, counter, label=paste("Knocking down parameter:", toString(counter), "/", toString(length(p_KD))))

    reacID_toKO <- which(grepl(model$namesSpecies[counter],model$reacID,fixed=TRUE))

    # model_KD$reacID <- model$reacID[-counter]
    model_KD$reacID <- model$reacID[-reacID_toKO]

    if (length(model_KD$reacID) > 0) {

      # Target_KOed <- strsplit(x = model$reacID[counter],split = "=")[[1]][2]
      Target_KOed <- model$namesSpecies[counter]

      if (model$reacID[1] == "EGF=PI3K") { # if CNOToy model
        estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=TRUE,HardConstraint=FALSE,Force=FALSE)
      } else if (model$reacID[1] == "PDGFL=PDGFR") { # if PDGF modeel
        # estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE,ORlist=c("Grb2SOS_OR_GabSOS=GGSOS"))
        estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE)
      } else {
        estim <- CNORprob_buildModel(CNOlist,model_KD,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE)
      }

      # If the KOed interaction left with target being not activated anymore -> Fix target to zero (instead of random number)

      if (sum(unique(nchar(estim$Interactions[,3]))==1)==1 && length(unique(nchar(estim$Interactions[,3])))==1) {
        if (sum(grepl(paste(Target_KOed,sep=""),estim$Interactions[,3],fixed = TRUE))==0) { # Somehow the grepl script doens't work properly for a single letter
          Target_KOed_Idx <- which(Target_KOed==estim$state_names)
          # if (sum(grepl(Target_KOed_Idx,estim$Input_vector[1,],fixed=TRUE))==0) { # If the KOed node is NOT an input
            estim$Input_vector <- cbind(estim$Input_vector,rep(0,dim(estim$Input_vector)[1]))
            estim$Input_index  <- cbind(estim$Input_index,rep(Target_KOed_Idx,dim(estim$Input_index)[1]))
          # }
        }
      } else {
        if (sum(grepl(paste("\\b",Target_KOed,"\\b",sep=""),estim$Interactions[,3],fixed = TRUE))==0) {
          Target_KOed_Idx <- which(Target_KOed==estim$state_names)
          if (sum(grepl(Target_KOed_Idx,estim$Input_index[1,],fixed=TRUE))==0) { # If the KOed node is NOT an input
            estim$Input_vector <- cbind(estim$Input_vector,rep(0,dim(estim$Input_vector)[1]))
            estim$Input_index  <- cbind(estim$Input_index,rep(Target_KOed_Idx,dim(estim$Input_index)[1]))
          }
            # Setting all downstream nodes of KOed node to be zero

            Reac2Remove <- model$reacID[reacID_toKO]
            AllTargets <- NULL
            for (counter2 in 1:length(Reac2Remove)) {
              AllTargets <- c(AllTargets,strsplit(Reac2Remove[counter2],split = "=")[[1]][2])
            }
            AllTargets <- unique(AllTargets)
            AllTargets_Idx <- NULL
            for (counter3 in 1:length(AllTargets)) {
              AllTargets_Idx <- c(AllTargets_Idx,which(AllTargets[counter3]==estim$state_names))
            }
            for (counter4 in 1:length(AllTargets)) {
              estim$Input_vector <- cbind(estim$Input_vector,rep(0,dim(estim$Input_vector)[1]))
              estim$Input_index  <- cbind(estim$Input_index,rep(AllTargets_Idx[counter4],dim(estim$Input_index)[1]))
            }

        }
      }
      estim$maxtime <- estim_orig$maxtime

      res <- NULL

      tryCatch(
        res <- CNORprob_optimise(estim,optRound=optRound_KO)
      ,
      error=function(e) print("Solver failed... continue the next optimisation round"))

      if (!is.null(res)) {

        CNORprob_plotFit(model_KD,CNOlist,estim,res,plotPDF=FALSE)

        SSE_L1 <- res$OptResults[1,1] # fetch the best cost from optimizations
        SSE <- SSE_L1 - (L1Reg*sum(res$ParamResults[,2]!=0))
        cost_KD[counter] <- SSE # fetch the best cost from optimizations

        # reduced model (-1 parameter)

        N_r = length(estim$Output_vector) # number of datapoints
        p_r= length(estim$param_vector) # number of parameters

        N_r_All <- c(N_r_All,N_r)
        p_r_All <- c(p_r_All,p_r)

        AIC_KD[counter] = N_r*log(cost_KD[counter]/N_r) + 2*p_r

      } else { # if the optimiser failed
        cost_KD[counter] <- NA
        AIC_KD[counter] <- NA
      }

    } else {
      print(paste("All interactions were removed from knocking down: Node ",model$namesSpecies[counter],sep=""))
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
  barplot(All_AIC, main="Node-KO analysis", names.arg = c("Orig",model$namesSpecies), col = cols, ylab = "AIC",cex.names = 0.8,las=2)
  abline(h=All_AIC[1],col="red")

  All_Cost <- c(SSE_orig,cost_KD)
  # All_KO_ParamNames <- c("Orig",estim$param_vector)
  All_KO_ParamNames <- c("Orig",model$reacID)
  cols <- c("gray","red")[(All_Cost < All_Cost[1])+1]
  cols[1] <- "blue" # First column always in blue
  barplot(All_Cost, main="Node-KO analysis", names.arg = c("Orig",model$namesSpecies), col = cols, ylab = "SSE",cex.names = 0.8,las=2)
  abline(h=All_Cost[1],col="red")

  # Exporting figures
  pdf("Results/Figures/CNORprob_nodeKO_AIC.pdf")
  barplot(All_AIC, main="Node-KO analysis", names.arg = c("Orig",model$namesSpecies), col = cols, ylab = "AIC",cex.names = 0.8,las=2)
  abline(h=All_AIC[1],col="red")
  dev.off()
  pdf("Results/Figures/CNORprob_nodeKO_Cost.pdf")
  barplot(All_Cost, main="Node-KO analysis", names.arg = c("Orig",model$namesSpecies), col = cols, ylab = "SSE",cex.names = 0.8,las=2)
  abline(h=All_Cost[1],col="red")
  dev.off()

  estim_Result$nodeKO$AIC_KD   <- All_AIC
  estim_Result$nodeKO$Cost_KD  <- All_Cost
  estim_Result                <<- estim_Result
  estim                       <<- estim_orig
  estim                       <<- estim

  return(estim_Result)

}

# ======= End of the script ======= #
