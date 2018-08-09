#' @title CNORprob_mapModel
#'
#' @description Map optimised parameters from CNORprob back to each reaction in model$reacID to be compatible with model plotting function in the CellNOptR package
#'
#' @export

CNORprob_mapModel = function(model,CNOlist,estim,res) {

  Splitted_reac <- list()
  Source_reac <- NULL
  Target_reac <- NULL
  for (counter in 1:length(model$reacID)) {
    Splitted_reac[[counter]] <- strsplit(model$reacID[counter],split = "=")
    Source_reac <- c(Source_reac,Splitted_reac[[counter]][[1]][1])
    Target_reac <- c(Target_reac,Splitted_reac[[counter]][[1]][2])
  }

  # BasalIdx <- which("basal"==Source_reac)
  #
  # options(warn = 1)
  #
  # if (length(BasalIdx)>0) {
  #   Source_reac <- Source_reac[-BasalIdx]
  #   Target_reac <- Target_reac[-BasalIdx]
  #
  #   model$reacID <- model$reacID[-BasalIdx]
  #   model$namesSpecies <- model$namesSpecies[-which("basal"==model$namesSpecies)]
  #   model$interMat <- model$interMat[,-BasalIdx]
  #   model$interMat <- model$interMat[-which("basal"==rownames(model$interMat)),]
  #   model$notMat <- model$notMat[,-BasalIdx]
  #   model$notMat <- model$notMat[-which("basal"==rownames(model$notMat)),]
  #
  #   CNOlist@cues <- as.matrix(CNOlist@cues[,-which("basal"==colnames(CNOlist@cues))])
  #   CNOlist@stimuli <- as.matrix(CNOlist@stimuli[,-which("basal"==colnames(CNOlist@stimuli))])
  #
  #   print("==================")
  #
  #   print("Note - The node 'basal' and all associated interactions were removed from visualisation")
  #
  #   print("==================")
  #
  # }

  print("==================")
  # warning("Note - Two separated interactions coming to a single node represent a 'competitive interaction'
  #         in CNORprob and not the 'OR' gate like the standard representation. Please refer to
  #         the original article of 'FALCON' which CNORprob was based on for more information.")

  print("Note - Two separated interactions coming to a single node represent a 'competitive interaction'")
  print("in CNORprob and not the 'OR' gate like the standard representation. Please refer to ")
  print("the original article of 'FALCON' which CNORprob was based on for more information.")


  print("==================")

  options(warn = 0)

  Unique_Targets <- unique(Target_reac)

  # Multiple entries statistic (for OR gate expansion)
  Target_reac_tab <- table(Target_reac)
  Multi_inputs_idx <- which(Target_reac_tab > 1)
  Multi_inputs_names <- names(Target_reac_tab[Multi_inputs_idx])
  Multi_inputs_freq <- as.numeric(Target_reac_tab[Multi_inputs_idx])

  OR_Reac_ToAdd <- NULL
  for (counter in 1:length(Multi_inputs_names)) {
    if (Multi_inputs_freq[counter]==2) {
      Current_Multi_inputs_names <- Source_reac[which(Multi_inputs_names[counter]==Target_reac)]
      if (sum(grepl("!",Current_Multi_inputs_names,fixed = TRUE))==0) {
        OR_Reac_ToAdd <- c(OR_Reac_ToAdd,paste(Current_Multi_inputs_names[1],"|",Current_Multi_inputs_names[2],"=",Multi_inputs_names[counter],sep=""))
      }
    }
  }

  bString <- NULL
  AND_NOT_Idx <- NULL
  AND_NOT_IntAct <- NULL
  AND_NOT_Names <- NULL

  for (counter in 1:length(model$reacID)) {


    # Clarify AND gate first
    if (grepl("+",Source_reac[counter],fixed = TRUE)) {
      # Debugged 08.07.18
      # Source_AND_Names_plus <- strsplit(Source_reac[counter],split = "+")
      # Source_AND_Names <- Source_AND_Names_plus[[1]][c(1,3)]
      Source_AND_Names <- strsplit(Source_reac[counter],split = "+",fixed=TRUE)
      Source_Idx <- NULL

      if (sum(grepl(pattern = "!",x = Source_AND_Names[[1]],fixed = T))==0) {
        grepl("+",Source_reac[counter],fixed = TRUE)

        for (counter_AND in 1:length(Source_AND_Names[[1]])) {
          Source_Idx <- c(Source_Idx,which(Source_AND_Names[[1]][counter_AND]==estim$Interactions[,1] & Target_reac[counter]==estim$Interactions[,3] & estim$Interactions[,5]=="A"))
        }
        IntAct_Idx <- Source_Idx[1] # Take either the 1st or 2nd index (should be the same)

      } else {

        AND_NOT_Idx <- c(AND_NOT_Idx, counter)
        AND_NOT_Names <- rbind(AND_NOT_Names, c(Source_AND_Names[[1]],Target_reac[counter]))
        # remove ! from the list of inputs (if any)
        for (counter_src in 1:length(Source_AND_Names[[1]])) {
          if (grepl(pattern = "!",x = Source_AND_Names[[1]][counter_src],fixed = T)) {
            Source_AND_Names[[1]][counter_src] <- substr(x = Source_AND_Names[[1]][counter_src],start = 2,stop = nchar(Source_AND_Names[[1]][counter_src]))
          }
        }
        for (counter_AND in 1:length(Source_AND_Names[[1]])) {
          Source_Idx <- c(Source_Idx,which((Source_AND_Names[[1]][counter_AND]==estim$Interactions[,1]) & Target_reac[counter]==estim$Interactions[,3]))
        }
        IntAct_Idx <- NA # Take either the 1st or 2nd index (should be the same)
        AND_NOT_IntAct <- rbind(AND_NOT_IntAct, Source_Idx)

      }


    } else {

      if (grepl("!",Source_reac[counter],fixed = TRUE)) {
        Source_Idx <- which(substring(Source_reac[counter],2)==estim$Interactions[,1])
      } else {
        Source_Idx <- which(Source_reac[counter]==estim$Interactions[,1])
      }
      Target_Idx <- which(Target_reac[counter]==estim$Interactions[,3])
      IntAct_Idx <- intersect(Source_Idx,Target_Idx)

      # If the different types of interaction are assigned to the same reaction
      if (length(IntAct_Idx)>1) {
        if(length(unique(estim$Interactions$V2[IntAct_Idx]))>1) {

          Sign_ReacID <- grepl("!",strsplit(model$reacID[counter],split = "="),fixed=TRUE) # get the sign
          if (Sign_ReacID==FALSE) { # if activation
            Positive_Idx <- which("->"==estim$Interactions[,2])
            IntAct_Idx <- intersect(IntAct_Idx,Positive_Idx)
          } else if (Sign_ReacID==TRUE) { # if inhibition
            Negative_Idx <- which("-|"==estim$Interactions[,2])
            IntAct_Idx <- intersect(IntAct_Idx,Negative_Idx)
          }

        } else {

          IdxToKeep <- which(estim$Interactions$V5[IntAct_Idx]=="N")
          IntAct_Idx <- IntAct_Idx[IdxToKeep]

        }
      }

      # Now check if there is any expanded OR interaction
      OR_exist <- sum(grepl("_OR_",estim$Interactions[Target_Idx,1],fixed = TRUE))>0
      if (OR_exist) {
        OR_current_params <- estim$Interactions[Target_Idx,4]
        OR_current_params_kn <- OR_current_params[grepl("kn",OR_current_params,fixed = TRUE)]
        OR_current_params_kn_map <- NULL
        for (counter2 in 1:length(OR_current_params_kn)) {
          OR_current_params_kn_map <- c(OR_current_params_kn_map,which(OR_current_params_kn[counter2]==estim$Interactions[,4]))
        }
        OR_current_params_kn_idx <- NULL
        for (counter3 in 1:length(OR_current_params_kn)) {
          OR_current_params_kn_idx <- c(OR_current_params_kn_idx,which(OR_current_params_kn[counter3]==estim$param_vector))
        }

        OR_current_params_kn_val <- res$BestFitParams[OR_current_params_kn_idx]

        # if OR param value is more than the others
        if (OR_current_params_kn_val[length(OR_current_params_kn_val)]==max(OR_current_params_kn_val)) {
          IntAct_Idx <- tail(OR_current_params_kn_map,n=1) # Then take the OR value (last one)
        }
      }
    }

    Current_Param <- estim$Interactions[IntAct_Idx,4]
    if (sum(!is.na(Current_Param))) {
      if (is.na(as.numeric(Current_Param)>0)) {
        Current_Param_Idx <- which(Current_Param==estim$param_vector)
        Current_Param_Val <- res$BestFitParams[Current_Param_Idx]
      } else {
        Current_Param_Val <- as.numeric(Current_Param)
      }
    } else {
      Current_Param_Val = NA
    }


    # print(paste0(counter," -- ",Current_Param_Val)) # Might be useful for debugging

    bString <- c(bString, Current_Param_Val)

  }


  # Map back all results from AND NOT gate (if any)
  if (!is.null(AND_NOT_Idx)) {
    for (counter_retake in 1:length(AND_NOT_Idx)) {
      for (counter_mapANDNOT in 1:2) {
        model$reacID <- c(model$reacID,paste0(AND_NOT_Names[counter_retake,counter_mapANDNOT],"=",AND_NOT_Names[counter_retake,3]))
        Current_Param <- estim$Interactions[AND_NOT_IntAct[counter_retake,counter_mapANDNOT],4]
        Current_Param_Idx <- which(Current_Param==estim$param_vector)
        Current_Param_Val <- res$BestFitParams[Current_Param_Idx]
        bString <- c(bString, Current_Param_Val)
      }
    }
    model$reacID <- model$reacID[-AND_NOT_Idx]
    bString <- bString[-AND_NOT_Idx]
  }

  MappedProb <- list(bString=bString,model=model,CNOlist=CNOlist)

  return(MappedProb)

}

# ======= End of the script ====== #
