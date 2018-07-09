#' @title CNORprob_simulate
#'
#' @description Simulate CNORprob model with a defined set of parameters: either 'initial guess' or 'optimised'
#'
#' @export

CNORprob_simulate = function(estim,SimParams,SimIterations) {

  # Plug in the best fit parameters for model simulation
  # k <- BestFitParams
  k <- SimParams
  PlotIterations <- SimIterations

  # Keep current parameter value for troubleshooting
  # estim$current_k     <<- k

  # Retrieve back all necessary variables from estim (global)
  estim$Interactions   -> Interactions

  # Need to turn Interaction back as a list
  Interactions_List <- list()
  for (counter in 1:nrow(Interactions)) {
    Interactions_List[[counter]] <- Interactions[counter,]
  }
  Interactions <- Interactions_List

  estim$SSthresh       -> SSthresh
  estim$Input_vector   -> Input_vector
  estim$Input_index    -> Input_index
  estim$Output_vector  -> Output_vector
  estim$Output_index   -> Output_index
  estim$SD_vector      -> SD_vector
  estim$ma             -> ma
  estim$mi             -> mi
  estim$state_names    -> state_names
  estim$param_vector   -> param_vector
  estim$param_index    -> param_index
  estim$NrStates       -> NrStates
  estim$NrParams       -> NrParams
  estim$FixBool        -> FixBool

  # Initialise collecting variables
  MeanStateValueAll <- matrix(NA,dim(estim$Output_vector)[1],length(state_names))
  StdStateValueAll  <- matrix(NA,dim(estim$Output_vector)[1],length(state_names))
  MeanCostAll       <- rep(NA,dim(estim$Output_vector)[1])
  StdCostAll        <- rep(NA,dim(estim$Output_vector)[1])

  n = length(state_names)

  if (n<=25) {
    initial_t=10
    step_t=10
  } else if (n>25 & n<=100) {
    initial_t=100
    step_t=10
  } else if (n>100) {
    initial_t=300
    step_t=150
  }

  # Assign sample params to ma and mi
  # HyperFALCON k-mapping

  process_param_index <- param_index
  process_FixBool <- FixBool

  All_params <- c()
  for (i in 1:length(Interactions)) {
    Extracted_params <- Interactions[[i]][4]
    All_params <- c(All_params,Extracted_params)
  }

  IdxFixVal  <- which(!is.na(as.numeric(All_params))) # Get indicies of fixed values
  if (length(IdxFixVal)>0) {
    process_param_index <- process_param_index[-c(IdxFixVal),] # Remove interactions with fixed values
    IdxFixVal_BL <- rep(FALSE,length(FixBool))
    IdxFixVal_BL[IdxFixVal] <- TRUE
    # if (intersect(IdxFixVal_BL,FixBool)) {
    if (sum(IdxFixVal_BL & FixBool)>0) { # Check if Fix Boolean value are picked
      process_FixBool <- process_FixBool[-c(IdxFixVal)] # Then remove Boolean interactions with fixed value
    }
  }

  # Re-evaluate Boolmax after constant parameter value removal
  if (!is.null(dim(process_param_index))) { # if not a single vector
    if (max(process_param_index[,5])>0) { # if there exist Boolmax in re-initialiseBool_Remained <- sum(process_param_index[,5] != 0)/2
      Bool_Remained <- sum(process_param_index[,5] != 0)/2
      Bool_Remained_Idx <- which(process_param_index[,5] != 0)
      Bool_Counter <- 1
      Bool_NewCount <- 1
      while (Bool_Counter <= Bool_Remained*2) {
        process_param_index[Bool_Remained_Idx[Bool_Counter],5] <- Bool_NewCount
        process_param_index[Bool_Remained_Idx[Bool_Counter+1],5] <- Bool_NewCount
        Bool_Counter <- Bool_Counter + 2
        Bool_NewCount <- Bool_NewCount + 1
      }
    }
  }

  if (is.null(dim(process_param_index))) {
    BoolMax=max(process_param_index[5])
  } else {
    BoolMax=max(process_param_index[,5])
  }

  if (BoolMax>0) {
    for (counter in 1:BoolMax) {
      BoolIdx <- which(process_param_index[,5]==counter)
      if (sum(process_FixBool[BoolIdx])>0) { # If Boolean gate is already fixed
        process_param_index <- process_param_index[-c(BoolIdx),] # Remove both
      } else { # If not
        process_param_index <- process_param_index[-c(BoolIdx[2]),] # Remove just the second entry
      }
    }
  }

  pd <- process_param_index

  if (is.null(dim(process_param_index))) {
    kInd <- 1
  } else {
    kInd <- 1:dim(process_param_index)[1]
  }

  kmap <- k[kInd] # extend k including boolean gates. Now k2 has the same length as param_index
  kmap <- kmap[!is.na(kmap)] # remove NA (if any)

  if (is.null(dim(process_param_index))) { # if param_index is a vector
    kA <- kmap[which(pd[3] == 1)]
    kI <- kmap[which(pd[4] == 1)]
    if (length(kA) > 0) {ma[pd[1],pd[2]]=kA}
    if (length(kI) > 0) {mi[pd[1],pd[2]]=kI}
  } else { # if param_index is a matrix
    kA <- kmap[which(pd[,3] == 1)]
    kI <- kmap[which(pd[,4] == 1)]
    for (i in 1:length(kA)) {ma[pd[which(pd[,3] == 1)[i],1],pd[which(pd[,3] == 1)[i],2]]=kA[i]}
    for (i in 1:length(kI)) {mi[pd[which(pd[,4] == 1)[i],1],pd[which(pd[,4] == 1)[i],2]]=kI[i]}
  }

  # Get Boolean gate indices in for-loop
  if (is.null(dim(param_index))) {
    BoolMax=max(param_index[5])
  } else {
    BoolMax=max(param_index[,5])
  }

  if (BoolMax>0) {

    Gate_fill_indices <- NULL
    Gate_value_fill <- NULL

    for (counter2 in 1:BoolMax) {
      gate_indices <- which(param_index[,5]==counter2)
      current_ma_value <- unique(ma[unique(param_index[gate_indices,1]),param_index[gate_indices,2]])
      Current_fill_indices <- c(unique(param_index[gate_indices,1]), param_index[gate_indices[1],2], param_index[gate_indices[2],2], unique(param_index[gate_indices,6]))
      Gate_value_fill=rbind(Gate_value_fill, current_ma_value)
      Gate_fill_indices=rbind(Gate_fill_indices, Current_fill_indices)
    }
    Temp_BoolVal=matrix(data = 0,nrow = n,ncol = dim(Output_index)[1])
  }

  Inputs = Input_vector
  Measurements = Output_vector

  diff <- 0

  for (counter_exp in 1:dim(Measurements)[1]) {

    Collect_x <- matrix(NA,PlotIterations,length(state_names))
    diff_ALL  <- rep(NA,PlotIterations)

    x <- runif(n,min=0,max=1)
    dim(x) <- c(length(x),1)
    x[Input_index[counter_exp,]] <- Inputs[counter_exp,]
    x_orig <- x

    xmeas <- Measurements[counter_exp,]
    mask <- is.nan(xmeas)
    if (sum(mask) > 0) {
      xmeas[mask] <- 0
    }
    for (counter_plot in 1:PlotIterations) {
      diff <- 0
      break_point_ss <- TRUE
      runs <- initial_t
      x <- x_orig

      while (break_point_ss==TRUE) {

        pre_x <- x

        for (counter in 1:runs) {

          if (BoolMax>0) {
            if (is.matrix(Gate_fill_indices)) {
              for (counter2 in 1:BoolMax) {
                if (Gate_fill_indices[counter2,4]==1) { # OR gate
                  Temp_BoolVal[counter2]=ma[Gate_fill_indices[counter2,1],Gate_fill_indices[counter2,2]] * (1-(1-(x[Gate_fill_indices[counter2,2]]))*(1-x[Gate_fill_indices[counter2,3]]));
                } else if (Gate_fill_indices[counter2,4]==2) { # AND gate
                  Temp_BoolVal[counter2]=ma[Gate_fill_indices[counter2,1],Gate_fill_indices[counter2,2]] * (x[Gate_fill_indices[counter2,2]]*x[Gate_fill_indices[counter2,3]]);
                }
              }
            } else { # Only one gate
              if (Gate_fill_indices[4]==1) { # OR gate
                Temp_BoolVal=ma[Gate_fill_indices[1],Gate_fill_indices[2]] * (1-(1-(x[Gate_fill_indices[2]]))*(1-x[Gate_fill_indices[3]]));
              } else if (Gate_fill_indices[4]==2) { # AND gate
                Temp_BoolVal=ma[Gate_fill_indices[1],Gate_fill_indices[2]] * (x[Gate_fill_indices[2]]*x[Gate_fill_indices[3]]);
              }
            }
          }

          x <- (ma %*% x) * (matrix(1,n,1)-(mi %*% x))

          if (BoolMax>0) {
            if (is.matrix(Gate_fill_indices)) {
              for (counter3 in 1:BoolMax) {
                x[Gate_fill_indices[counter3,1]]=Temp_BoolVal[counter3]
              }
            } else { # Only one gate
              x[Gate_fill_indices[1]]=Temp_BoolVal
            }
          }

        }

        if (sum(abs(pre_x-x))>SSthresh) {
          runs = runs+step_t
        } else {
          break_point_ss=FALSE
        }

      }
      # mask <- is.nan(xmeas)
      if (sum(mask) > 0) {
        xmeas[mask] <- x[Output_index[counter_exp,]][mask]
        # x[Output_index[counter_exp,]][mask] <- 0
        # xmeas[mask] <- 0
      }
    diff=diff+(sum((x[Output_index[counter_exp,]]-xmeas)^2/length(estim$Output_vector)))
    Collect_x[counter_plot,] <- x
    diff_ALL[counter_plot] <- diff
    }


  MeanAllState  <- colMeans(Collect_x)
  StdAllState   <- apply(Collect_x,2,sd)
  mean_diff_ALL <- mean(diff_ALL)
  std_diff_ALL  <- sd(diff_ALL)

  MeanStateValueAll[counter_exp,] <- MeanAllState
  StdStateValueAll[counter_exp,]  <- StdAllState
  MeanCostAll[counter_exp]        <- mean_diff_ALL
  StdCostAll[counter_exp]         <- std_diff_ALL

  }

  SimResAll <- list()
  SimResAll$MeanStateValueAll <- MeanStateValueAll
  SimResAll$StdStateValueAll  <- StdStateValueAll
  SimResAll$MeanCostAll       <- MeanCostAll
  SimResAll$StdCostAll        <- StdCostAll

  return(SimResAll)

}

# ======= End of the script ====== #
