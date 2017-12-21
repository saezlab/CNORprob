#' @title CNORprob_writeConstraint
#'
#' @description Write constraints for probabilistic logic models before passing to the optimiser
#'
#' @export

CNORprob_writeConstraint = function(Interactions,HLbound,state_names,HardConstraint,LPSA=NULL) {

  All_states <- c()
  colnames(Interactions) <- NULL
  rownames(Interactions) <- NULL
  Interactions_Matrix <- Interactions
  Interactions_List <- list()
  for (counter in 1:nrow(Interactions)) {
    Interactions_List[[counter]] <- Interactions[counter,]
  }
  Interactions <- Interactions_List

  # for (i in 1:length(Interactions)) {
  #   Extracted_states_inputs <- c(Interactions[[i]][1])
  #   Extracted_states <- c(Extracted_states_inputs)
  #   All_states <- c(All_states,Extracted_states)
  # }
  # for (i in 1:length(Interactions)) {
  #   Extracted_states_outputs <- c(Interactions[[i]][3])
  #   Extracted_states <- c(Extracted_states_outputs)
  #   All_states <- c(All_states,Extracted_states)
  # }
  # state_names <- unique(All_states)

  ma=matrix(0,length(state_names),length(state_names))
  mi=matrix(0,length(state_names),length(state_names))

  All_params <- c()
  for (i in 1:length(Interactions)) {
    Extracted_params <- Interactions[[i]][4]
    All_params <- c(All_params,Extracted_params)
  }
  param_names <- unique(All_params)

  # Remove fixed parameter value from the list of parameters
  break_point=TRUE;
  count=1;
  while (break_point) {
    if (count<=length(param_names)) {
      if (!is.na(as.numeric(param_names[count]))) {
        # param_names=param_names[param_names!=param_names[count]];
        param_names=param_names[-count];
      } else {
        count=count+1;
      }
    } else {
      break_point=0;
    }
  }
  param_vector<-param_names;
  dim(param_vector) <- c(length(param_vector),1)

  # 2nd step: Read in interaction and define parameters in ma
  param_index=c()
  BoolSet=1
  RemovedIntAct <- NULL

  for (counter in 1:length(Interactions)) {

    idx_input=match(Interactions[[counter]][1],state_names)
    idx_output=match(Interactions[[counter]][3],state_names)

    if (length(Interactions[[counter]])==5) {
      HLcon <- NULL
    } else if (length(Interactions[[counter]])==6) {
      if (grepl(Interactions[[counter]][6],"D",fixed=TRUE)) {
        HLcon <- 0;
      } else if (grepl(Interactions[[counter]][6],"H",fixed=TRUE)) {
        HLcon <- 1;
      } else if (grepl(Interactions[[counter]][6],"L",fixed=TRUE)) {
        HLcon <- -1;
      }
    }

    if (grepl(Interactions[[counter]][2],"->",fixed=TRUE)) { # if activates
      ma[idx_output,idx_input]=1

      if (grepl(Interactions[[counter]][5],"O",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 1, 0, BoolSet, 1, HLcon))
      } else if (grepl(Interactions[[counter]][5],"A",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 1, 0, BoolSet, 2, HLcon))
      } else if (grepl(Interactions[[counter]][5],"N",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 1, 0, 0, 0, HLcon))
      }

      # Map assigned (real) parameter value & remove entry from param_index
      if (!is.na(as.numeric(Interactions[[counter]][4])>0)) {
        ma[idx_output,idx_input]=as.numeric(Interactions[[counter]][4])
        if (as.character(Interactions[[counter]][5])=="N") {
          # param_index <- param_index[-c(nrow(param_index)),]
          RemovedIntAct <- c(RemovedIntAct,counter)
        }
      }

    } else if (grepl(Interactions[[counter]][2],"-|",fixed=TRUE)) { # if inhibits
      mi[idx_output,idx_input]=1;

      if (grepl(Interactions[[counter]][5],"O",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 0, 1, BoolSet, 1, HLcon))
      } else if (grepl(Interactions[[counter]][5],"A",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 0, 1, BoolSet, 2, HLcon))
      } else if (grepl(Interactions[[counter]][5],"N",fixed=TRUE)) {
        param_index=rbind(param_index, c(idx_output, idx_input, 0, 1, 0, 0, HLcon))
      }

      # Map assigned (real) parameter value & remove entry from param_index
      if (!is.na(as.numeric(Interactions[[counter]][4])>0)) {
        mi[idx_output,idx_input]=as.numeric(Interactions[[counter]][4])
        if (as.character(Interactions[[counter]][5])=="N") {
          # param_index <- param_index[-c(nrow(param_index)),]
          RemovedIntAct <- c(RemovedIntAct,counter)
        }
      }
    }
  }
  rownames(param_index) <- NULL

  # Assign Boolean indicator

  unique_output=unique(param_index[,1])
  BoolSet=1;

  for (counter in 1:length(unique_output)) {
    BoolSet_index <- which(param_index[,1]==unique_output[counter])
    if (length(BoolSet_index)>1) {
      BoolSet_act_index <- which(param_index[BoolSet_index,3]==1 & param_index[BoolSet_index,5]==1)
      BoolSet_inh_index <- which(param_index[BoolSet_index,4]==1 & param_index[BoolSet_index,5]==1)
      # BoolSet_act_index <- which(param_index[BoolSet_index,3]==1)
      # BoolSet_inh_index <- which(param_index[BoolSet_index,4]==1)
      if ((length(BoolSet_act_index)>1) & (sum(unique(param_index[BoolSet_index[BoolSet_act_index],6])) != 0)) {
        param_index[BoolSet_index[BoolSet_act_index],5]=BoolSet;
        BoolSet=BoolSet+1;
      }
      if ((length(BoolSet_inh_index)>1) & (sum(unique(param_index[BoolSet_index[BoolSet_inh_index],6])) != 0)) {
        param_index[BoolSet_index[BoolSet_inh_index],5]=BoolSet;
        BoolSet=BoolSet+1;
      }
    }
  }

  # Identify fixed Boolean indices

  # FixVal  <- as.numeric(All_params)>0
  FixVal  <- !is.na(as.numeric(All_params))
  BoolPos <- param_index[,5]>0
  # FixBool <- !is.na(FixVal == BoolPos)
  FixBool <- FixVal & BoolPos

  # 3rd step: Find inputs indices (be in input column AND NOT in output column) & assign to ma

  All_sources <- c()
  for (i in 1:length(Interactions)) {
    Extracted_sources <- as.character(c(Interactions[[i]][1]))
    All_sources <- c(All_sources,Extracted_sources)
  }
  sources_names <- sort(unique(All_sources))

  All_targets <- c()
  for (i in 1:length(Interactions)) {
    Extracted_targets <- as.character(c(Interactions[[i]][3]))
    All_targets <- c(All_targets,Extracted_targets)
  }
  targets_names <- sort(unique(All_targets))

  source_idx=match(sources_names,state_names)
  target_idx=match(targets_names,state_names)

  Inputs_Indices <- source_idx[is.na(!match(source_idx,target_idx))]

  for (counter in 1:length(Inputs_Indices)) {
    ma[Inputs_Indices[counter],Inputs_Indices[counter]]=1;
  }

  # 4th step: Assign linear and non-linear constraints
  Aeq=NULL; beq=NULL; A=NULL; b=NULL;

  orig_param_index <- param_index # keep original param index
  orig_FixBool <- FixBool

  for (i in rev(1:length(Interactions))) {
    if (!is.na(as.numeric(Interactions[[i]][4])>0)) {
      param_index <- param_index[-c(i),]
      FixBool <- FixBool[-c(i)]
    }
  }

  # Re-evaluate Boolmax after constant parameter value removal
  if (!is.null(dim(param_index))) { # if not a single vector
    if (max(param_index[,5])>0) { # if there exist Boolmax in re-initialise
      Bool_Remained <- sum(param_index[,5] != 0)/2
      Bool_Remained_Idx <- which(param_index[,5] != 0)
      Bool_Counter <- 1
      Bool_NewCount <- 1
      while (Bool_Counter <= Bool_Remained*2) {
        param_index[Bool_Remained_Idx[Bool_Counter],5] <- Bool_NewCount
        param_index[Bool_Remained_Idx[Bool_Counter+1],5] <- Bool_NewCount
        Bool_Counter <- Bool_Counter + 2
        Bool_NewCount <- Bool_NewCount + 1
      }
    }
  }

  # remove duplicate entry with logic gates
  if (is.null(dim(param_index))) {
    BoolMax=max(param_index[5])
  } else {
    BoolMax=max(param_index[,5])
  }

  if (BoolMax>0) {
    for (counter in 1:BoolMax) {
      BoolIdx <- which(param_index[,5]==counter)
      if (sum(FixBool[BoolIdx])>0) { # If Boolean gate is already fixed
        param_index <- param_index[-c(BoolIdx),] # Remove both
      } else { # If not
        param_index <- param_index[-c(BoolIdx[2]),] # Remove just the second entry
      }
    }
  }

  # Assign constraints

  if (!is.null(dim(param_index))) {
    for (counter in 1:length(unique_output)) {
      current_index <- which(param_index[,1]==unique_output[counter])
      if (length(current_index)>1) {
        act_index <- which(param_index[current_index,3]==1)
        inh_index <- which(param_index[current_index,4]==1)
        if (length(act_index)>1) {
        # if (length(act_index)>1 & isTRUE(unique(param_index[current_index[act_index],6])==0)) {
          if (HardConstraint==TRUE) {
            current_ma <- current_index[act_index]
            Aeq_temp <- matrix(0,1,length(param_names))
            Aeq_temp[current_ma] <- 1
            Aeq <- rbind(Aeq, Aeq_temp)
            beq <- rbind(beq, 1)
          } else if (HardConstraint==FALSE) {
            current_ma <- current_index[act_index]
            A_temp <- matrix(0,1,length(param_names))
            A_temp[current_ma] <- 1
            A <- rbind(A, A_temp)
            b <- rbind(b, 1)
          }
        }
        if (length(inh_index)>1) {
        # if (length(inh_index)>1 & isTRUE(unique(param_index[current_index[inh_index],6])==0)) {
          current_mi=current_index[inh_index]
          A_temp <- matrix(0,1,length(param_names))
          A_temp[current_mi]=1
          A <- rbind(A, A_temp)
          b <- rbind(b, 1)
        }
      }
    }
  }

  # Derive additional High/Low constraints with bounds

  LB <- matrix(0,1,length(param_vector))
  UB <- matrix(1,1,length(param_vector))

  if (length(Interactions[[1]])==6) {

    if (is.null(dim(param_index))) {
      if (param_index[7]==0) {
        LB=0; UB=1;
      } else if (param_index[7]==1) {
        LB=HLbound; UB=1;
      } else if (param_index[7]==-1) {
        LB=0; UB=HLbound;
      }
    } else {
      for (counter in 1:dim(param_index)[1]) {
        if (param_index[counter,7]==0) {
          LB[counter]=0; UB[counter]=1;
        } else if (param_index[counter,7]==1) {
          LB[counter]=HLbound; UB[counter]=1;
        } else if (param_index[counter,7]==-1) {
          LB[counter]=0; UB[counter]=HLbound;
        }
      }
    }
  }

  param_index=orig_param_index; # return original param_index
  FixBool <- orig_FixBool; # return original FixBool index

  if (!is.null(LPSA)) {
    param_index <- param_index[-LPSA,]
    FixBool <- FixBool[-LPSA]
  }

  return(list(ma=ma,mi=mi,A=A,b=b,Aeq=Aeq,beq=beq,state_names=state_names,param_vector=param_vector,param_index=param_index,LB=LB,UB=UB,FixBool=FixBool))

}

# ======= End of the script ======= #

