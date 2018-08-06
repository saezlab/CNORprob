#' @title preprocessing_Prob
#'
#' @description preprocessing_Prob
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

preprocessing_Prob<-function(data=NULL, model, cutNONC=TRUE, compression=TRUE,
    expansion=FALSE, ignoreList=NA, maxInputsPerGate=2,verbose=TRUE){

    # why not doing this check here ? Does not cost too much
    if (is.null(data)!=TRUE){
        checkSignals(CNOlist=data,model=model)
    }

    # a copy of the model
    cutModel <- model

    if (cutNONC==TRUE && is.null(data)!=TRUE){
        # Find the indices, in the model, of the species that are inh/stim/sign
        indices<-indexFinder(CNOlist=data, model=model,    verbose=verbose)

        # Find the indices of the non-osb/non-contr
        temp_indices <- findNONC(model=model, indexes=indices,verbose=verbose)
        # Cut the nonc off the model
        cutModel <-cutNONC(model=model, NONCindexes=temp_indices)
    }

    if (compression == TRUE && is.null(data)!=TRUE){
        # Recompute the indices
        temp_indices<-indexFinder(CNOlist=data, model=cutModel)

        # Compress the model
        cutModel<-compressModel(model=cutModel,indexes=temp_indices)
    }

    # Recompute the indices. We can do it now because the expanson gate does not
    # remove species but only add and gates.
    #if (is.null(data)!=TRUE){
    #    indices<-indexFinder(CNOlist=data,model=cutModel)
    #}
    #else{
    #    indices <- NULL
    #}

    # Expand the gates
    if (expansion == TRUE){
        cutModel <- expandGates(model=cutModel, ignoreList=ignoreList,maxInputsPerGate=maxInputsPerGate)
    }


    # !!! CNORprob specific !!! #
    # 1. AND NOT gate in model$reacID will be automatically splitted into 2 separated interaction

    ToAdd_ANDreac <- NULL
    Splitted_ANDchar <- NULL
    Splitted_ANDcharEQ <- NULL
    ANDreac_Idx <- NULL
    for (counter in 1:length(cutModel$reacID)) {
      if (grepl("+",cutModel$reacID[counter],fixed = TRUE)) {
        # Splitted_ANDchar[[counter]] <- strsplit(model$reacID[counter],split = "+") # Fixed 05.07.18
        Splitted_ANDchar[[counter]] <- strsplit(cutModel$reacID[counter],split = "+",fixed = T)
        Splitted_ANDcharEQ[[counter]] <- strsplit(Splitted_ANDchar[[counter]][[1]][2],split = "=",fixed = T)

        # Check for AND NOT gate -> equal to default combined positive/negative gate in CNORprob
        if (length(grep("!",Splitted_ANDchar[[counter]][[1]][1],fixed=T))==0 & length(grep("!",Splitted_ANDcharEQ[[counter]][[1]][1],fixed=T))==0) {
          # ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste("&",Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep="")),
          #                        paste("&",Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep=""))
          # ANDreac_Idx <- c(ANDreac_Idx,counter)

        } else {
          if (length(grep("!",Splitted_ANDchar[[counter]][[1]][1],fixed=T))!=0) {
            ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste(Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep=""),
                                                       paste(Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep="")))
          } else {
            ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste(Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep=""),
                                                       paste(Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep="")))
          }
          print("======================================================")
          ANDreac_Idx <- c(ANDreac_Idx,counter)
          print(paste0("Note: AND NOT gate detected for: ", cutModel$reacID[counter]))
          print(paste0(cutModel$reacID[counter], " is now splitted into ",
                       ToAdd_ANDreac[length(ToAdd_ANDreac)]," and ",ToAdd_ANDreac[length(ToAdd_ANDreac)-1]))
          print("======================================================")

        }
      }
    }

    # Remove the original reacID with "+" gate and replace it with new splitted reactions
    if (!is.null(ToAdd_ANDreac)) {
      cutModel$reacID <- c(cutModel$reacID[-ANDreac_Idx],ToAdd_ANDreac)
    } else {
      cutModel$reacID <- cutModel$reacID
    }


    # !!! CNORprob specific !!! #
    # 2. If a node is only inhibited in the original model -> add basal node

    # BasalAdded=0
    # Splitted_reac <- list()
    # Source_reac <- NULL
    # Target_reac <- NULL
    # for (counter in 1:length(cutModel$reacID)) {
    #   Splitted_reac[[counter]] <- strsplit(cutModel$reacID[counter],split = "=")
    #   Source_reac <- c(Source_reac,Splitted_reac[[counter]][[1]][1])
    #   Target_reac <- c(Target_reac,Splitted_reac[[counter]][[1]][2])
    # }
    #
    # TargetTable <- table(Target_reac)
    # for (counter in 1:length(TargetTable)) {
    #   if (TargetTable[counter]==1) {
    #     if (length(grep("!", Source_reac[which(Target_reac==names(TargetTable[counter]))],fixed=T))>0) {
    #
    #       BasalAdded=BasalAdded+1
    #
    #       cutModel$reacID <- c(cutModel$reacID,paste0("basal=",names(TargetTable[counter])))
    #       cutModel$interMat <- cbind(cutModel$interMat,rep(0,nrow(cutModel$interMat)))
    #       colnames(cutModel$interMat)[ncol(cutModel$interMat)] <- paste0("basal=",names(TargetTable[counter]))
    #       cutModel$notMat <- cbind(cutModel$notMat,rep(0,nrow(cutModel$notMat)))
    #       colnames(cutModel$notMat)[ncol(cutModel$notMat)] <- paste0("basal=",names(TargetTable[counter]))
    #
    #       if (BasalAdded==1) { # Add basal only once for rownames
    #         cutModel$namesSpecies <- c(cutModel$namesSpecies,"basal")
    #         cutModel$interMat <- rbind(cutModel$interMat,rep(0,ncol(cutModel$interMat)))
    #         rownames(cutModel$interMat)[nrow(cutModel$interMat)] <- "basal"
    #         cutModel$notMat <- rbind(cutModel$notMat,rep(0,ncol(cutModel$notMat)))
    #         rownames(cutModel$notMat)[nrow(cutModel$notMat)] <- "basal"
    #       }
    #
    #       cutModel$interMat[which(rownames(cutModel$interMat)=="basal"),
    #                         which(colnames(cutModel$interMat)==paste0("basal=",names(TargetTable[counter])))] <- -1
    #       cutModel$interMat[which(rownames(cutModel$interMat)==names(TargetTable[counter])),
    #                         which(colnames(cutModel$interMat)==paste0("basal=",names(TargetTable[counter])))] <- 1
    #
    #
    #
    #       print("======================================================")
    #       print(paste0("Note: Single inhibited node detected for: ", names(TargetTable[counter])))
    #       print(paste0("The interaction 'basal=", names(TargetTable[counter]), "' has been added."))
    #
    #       print("======================================================")
    #
    #
    #     }
    #   }
    # }
    #
    # if (BasalAdded>0) {
    #   data@cues <- cbind(data@cues,rep(1,nrow(data@cues)))
    #   colnames(data@cues)[ncol(data@cues)] <- "basal"
    #   data@stimuli <- cbind(data@stimuli,rep(1,nrow(data@stimuli)))
    #   colnames(data@stimuli)[ncol(data@stimuli)] <- "basal"
    # }

    # since version 1.3.28 return only model, indices are recomputed in other
    # functions

    PreProcessResults <- list(cutModel=cutModel,data=data)

    return(PreProcessResults)
}
