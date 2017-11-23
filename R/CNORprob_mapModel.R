CNORprob_mapModel = function(model,estim,res) {

  Splitted_reac <- list()
  Source_reac <- NULL
  Target_reac <- NULL
  for (counter in 1:length(model$reacID)) {
    Splitted_reac[[counter]] <- strsplit(model$reacID[counter],split = "=")
    Source_reac <- c(Source_reac,Splitted_reac[[counter]][[1]][1])
    Target_reac <- c(Target_reac,Splitted_reac[[counter]][[1]][2])
  }
  
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
  for (counter in 1:length(model$reacID)) {
    
    # Clarify AND gate first
    if (grepl("+",Source_reac[counter],fixed = TRUE)) {
      Source_AND_Names_plus <- strsplit(Source_reac[counter],split = "+")
      Source_AND_Names <- Source_AND_Names_plus[[1]][c(1,3)]
      Source_Idx <- NULL
      for (counter_AND in 1:length(Source_AND_Names)) {
        Source_Idx <- c(Source_Idx,which(Source_AND_Names[counter_AND]==estim$Interactions[,1]))
      }
      IntAct_Idx <- Source_Idx[1] # Take either the 1st or 2nd index (should be the same)
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
        Sign_ReacID <- grepl("!",strsplit(model$reacID[counter],split = "="),fixed=TRUE) # get the sign
        if (Sign_ReacID==FALSE) { # if activation
          Positive_Idx <- which("->"==estim$Interactions[,2])
          IntAct_Idx <- intersect(IntAct_Idx,Positive_Idx)
        } else if (Sign_ReacID==TRUE) { # if inhibition
          Negative_Idx <- which("-|"==estim$Interactions[,2])
          IntAct_Idx <- intersect(IntAct_Idx,Negative_Idx)
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
    if (is.na(as.numeric(Current_Param)>0)) {
      Current_Param_Idx <- which(Current_Param==estim$param_vector)
      Current_Param_Val <- res$BestFitParams[Current_Param_Idx]
    } else {
      Current_Param_Val <- as.numeric(Current_Param)
    }

    bString <- c(bString, Current_Param_Val)
  }
  
  return(bString)
 
}

# ======= End of the script ====== #
 