#' @title CNORprob_buildModel
#'
#' @description Build a probablistic logic model description from pknmodel and CNOlist
#'
#' @export

CNORprob_buildModel = function(CNOlist,model,expandOR=FALSE,HardConstraint=TRUE,Force=TRUE,ORlist=NULL,L1Reg=0.01,HLbound=0.5,SSthresh=2e-16,
                               PlotIterations=1,rsolnp_options=list(rho=1,outer.iter=400,inner.iter=800,delta=1e-7,tol=1e-8,trace=1)) {

  # ===== Step 1 : Extract data from CNOlist  ==== #

  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }

  NrExp <- dim(getCues(CNOlist))[1]
  Inhibitor_names <- colnames(getInhibitors(CNOlist))

  if (is.null(Inhibitor_names)) {
    state_names <- model$namesSpecies
  } else {
    state_names <- c(model$namesSpecies,paste(Inhibitor_names,"i",sep=""))
    Inhibitor_names_with_i <- paste(Inhibitor_names,"i",sep="")
  }

  Stimuli_names <- colnames(getStimuli(CNOlist))
  if (is.null(Stimuli_names)) {
    stop(print('At least one stimuli is needed in MIDAS file'))
  }

  # if (is.null(Inhibitor_names)) {
  #   Input_names <- Stimuli_names
  # } else {
  #   Input_names <- c(Stimuli_names,Inhibitor_names_with_i)
  # }

  Input_vector <- getCues(CNOlist)
  # colnames(Input_vector) <- Input_names

  if (!is.null(Inhibitor_names)) {
    for (counter in 1:length(Inhibitor_names)) {
      colnames(Input_vector)[grep(Inhibitor_names[counter],colnames(Input_vector),fixed = T)] <- paste0(colnames(Input_vector)[grep(Inhibitor_names[counter],colnames(Input_vector),fixed = T)],"i")
    }
  }

  Input_names <- colnames(Input_vector)

  Input_index <- NULL
  for (counter in 1:length(Input_names)) {
    Input_index <- c(Input_index,which(Input_names[counter]==state_names))
  }
  Input_index <- rep.row(Input_index,dim(getCues(CNOlist))[1])


  Output_names <- colnames(getSignals(CNOlist)[[length(getSignals(CNOlist))]])
  Output_index <- NULL
  for (counter in 1:length(Output_names)) {
    Output_index <- c(Output_index,which(Output_names[counter]==state_names))
  }
  Output_index <- rep.row(Output_index,dim(getCues(CNOlist))[1])
  Output_vector <- getSignals(CNOlist)[[length(getSignals(CNOlist))]]

  SD_vector <- getVariances(CNOlist)[[length(getVariances(CNOlist))]]

  # ===== Step 2 Extract interactions from model -> build "Interaction" matrix ==== #

  # Preprocess reacID to split AND gate and assign special symbol
  ToAdd_ANDreac <- NULL
  Splitted_ANDchar <- NULL
  Splitted_ANDcharEQ <- NULL
  ANDreac_Idx <- NULL
  for (counter in 1:length(model$reacID)) {
    if (grepl("+",model$reacID[counter],fixed = TRUE)) {
      # Splitted_ANDchar[[counter]] <- strsplit(model$reacID[counter],split = "+") # Fixed 05.07.18
      Splitted_ANDchar[[counter]] <- strsplit(model$reacID[counter],split = "+",fixed = T)
      Splitted_ANDcharEQ[[counter]] <- strsplit(Splitted_ANDchar[[counter]][[1]][2],split = "=",fixed = T)
      # Check for AND NOT gate -> equal to default combined positive/negative gate in CNORprob
      if (length(grep("!",Splitted_ANDchar[[counter]][[1]][1],fixed=T))==0 & length(grep("!",Splitted_ANDcharEQ[[counter]][[1]][1],fixed=T))==0) {
        ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste("&",Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep="")),
                               paste("&",Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep=""))
        ANDreac_Idx <- c(ANDreac_Idx,counter)
      }
      # =========================================================== #

      # This part of the function has been moved to preprocessing_Prob.R instead

      # else {
      #   if (length(grep("!",Splitted_ANDchar[[counter]][[1]][1],fixed=T))!=0) {
      #     ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste(Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep=""),
      #                                                paste(Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep="")))
      #   } else {
      #     ToAdd_ANDreac <- rbind(ToAdd_ANDreac,rbind(paste(Splitted_ANDchar[[counter]][[1]][1],"=", Splitted_ANDcharEQ[[counter]][[1]][2],sep=""),
      #                                                paste(Splitted_ANDcharEQ[[counter]][[1]][1],"=",Splitted_ANDcharEQ[[counter]][[1]][2],sep="")))
      #   }
      #   ANDreac_Idx <- c(ANDreac_Idx,counter)
      # }
      #
      # =========================================================== #

    }
  }

  # Remove the original reacID with "+" gate and replace it with new splitted reactions
  if (!is.null(ToAdd_ANDreac)) {
    model$reacID_pl <- c(model$reacID[-ANDreac_Idx],ToAdd_ANDreac)
  } else {
    model$reacID_pl <- model$reacID
  }


  # Start preparing the internal interaction file for CNORprob
  Splitted_reac <- list()
  Source_reac <- NULL
  Target_reac <- NULL
  for (counter in 1:length(model$reacID_pl)) {
    Splitted_reac[[counter]] <- strsplit(model$reacID_pl[counter],split = "=")
    Source_reac <- c(Source_reac,Splitted_reac[[counter]][[1]][1])
    Target_reac <- c(Target_reac,Splitted_reac[[counter]][[1]][2])
  }
  Unique_Targets <- unique(Target_reac)

  Interactions <- as.data.frame(matrix(NA,length(model$reacID_pl)+length(Inhibitor_names),6))
  Interactions[,5] <- 'N'
  Interactions[,6] <- 'D'

  i <- 1
  j <- 0

  for (counter in 1:length(Unique_Targets)) {

    ReacIdx <- which(Target_reac == Unique_Targets[counter])

    if (expandOR) {

      if (length(ReacIdx) > 1) { # if more than 1 interaction per output
        ActIdx <- which(!grepl("!",Source_reac[ReacIdx],fixed = TRUE))
        InhIdx <- which(grepl("!",Source_reac[ReacIdx],fixed = TRUE))

        if (length(ActIdx) > 2 || length(InhIdx) > 2) {

          # Need to separate interactions into pairs here
          stop("The current CNORprob pipeline doesn't support multiple gates assignment.
          Please consider revising interaction files using a single logical gate e.g.
          [A 1 T;B 1 T;C 1 T] -> [A 1 AB;B 1 AB;AB 1 T;C 1 T]
          AND Please do not use expandGate during pre-processing step in CNORprob pipeline")

        }

        if (length(ActIdx) == 2) {

          for (counter2 in 1:length(ActIdx)) {
            if (grepl("&",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
              Interactions[i,1] <- substring(Source_reac[ReacIdx[counter2]],2)
            } else {
              Interactions[i,1] <- Source_reac[ReacIdx[counter2]]
            }
            Interactions[i,2] <- "->"
            Interactions[i,3] <- Target_reac[ReacIdx[counter2]]
            Interactions[i,4] <- paste('k',toString(counter+j),sep="")
            if (grepl("&",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
              Interactions[i,5] <- "A"
            } else {
              Interactions[i,5] <- "O"
            }
            i <- i+1
          }
        } else if (length(ActIdx) == 1) {
          j <- j+1
          Interactions[i,1] <- Source_reac[ReacIdx[ActIdx]]
          Interactions[i,2] <- "->"
          Interactions[i,3] <- Target_reac[ReacIdx[ActIdx]]
          Interactions[i,4] <- paste('k',toString(counter+j),sep="")
          i <- i+1
        }

        if (length(InhIdx) == 2) {

          for (counter2 in 1:length(ActIdx)) {
            if (grepl("&",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
              Interactions[i,1] <- substring(Source_reac[ReacIdx[counter2]],2)
            } else {
              Interactions[i,1] <- Source_reac[ReacIdx[counter2]]
            }
            Interactions[i,2] <- "-|"
            Interactions[i,3] <- Target_reac[ReacIdx[counter2]]
            Interactions[i,4] <- paste('k',toString(counter+j),sep="")
            if (grepl("&",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
              Interactions[i,5] <- "A"
            } else {
              Interactions[i,5] <- "O"
            }
            i <- i+1
          }
        } else if (length(InhIdx) == 1) {
          j <- j+1
          Interactions[i,1] <- substring(Source_reac[ReacIdx[InhIdx]],2)
          Interactions[i,2] <- "-|"
          Interactions[i,3] <- Target_reac[ReacIdx[InhIdx]]
          Interactions[i,4] <- paste('k',toString(counter+j),sep="")
          i <- i+1
        }

      } else {
        if (grepl("!",Source_reac[ReacIdx],fixed = TRUE)) {
          Interactions[i,1] <- substring(Source_reac[ReacIdx],2)
          Interactions[i,2] <- "-|"
        } else {
          Interactions[i,1] <- Source_reac[ReacIdx]
          Interactions[i,2] <- "->"
        }
        Interactions[i,3] <- Target_reac[ReacIdx]
        Interactions[i,4] <- paste('k',toString(counter+j),sep="")
        i <- i+1
      }
    } else {

      if (length(ReacIdx) > 1) { # if more than 1 interaction per output

        # ActIdx <- which(!grepl("!",Source_reac[ReacIdx],fixed = TRUE))
        # InhIdx <- which(grepl("!",Source_reac[ReacIdx],fixed = TRUE))
        # 
        # if (length(ActIdx) > 2 || length(InhIdx) > 2) {
        #   
        #   # Need to separate interactions into pairs here
        #   stop("The current CNORprob pipeline doesn't support multiple gates assignment.
        #        Please consider revising interaction files using a single logical gate e.g.
        #        [A 1 T;B 1 T;C 1 T] -> [A 1 AB;B 1 AB;AB 1 T;C 1 T]
        #        AND Please do not use expandGate during pre-processing step in CNORprob pipeline")
        #   
        # }

        for (counter2 in 1:length(ReacIdx)) {
          if (grepl("!",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
            Interactions[i,1] <- substring(Source_reac[ReacIdx[counter2]],2)
            Interactions[i,2] <- "-|"
          } else {
            Interactions[i,1] <- Source_reac[ReacIdx[counter2]]
            Interactions[i,2] <- "->"
          }
          Interactions[i,3] <- Target_reac[ReacIdx[counter2]]
          if (grepl("&",Source_reac[ReacIdx[counter2]],fixed = TRUE)) {
            Interactions[i,1] <- substring(Source_reac[ReacIdx[counter2]],2)
            Interactions[i,5] <- "A"
            if (counter2 != 1) { # only the first interaction get a new parameter name
              j <- j-1 # use the same index for all AND gate
            }
          }
          Interactions[i,4] <- paste('k',toString(counter+j),sep="")
          i <- i+1
          j <- j+1
        }

      } else {
        if (grepl("!",Source_reac[ReacIdx],fixed = TRUE)) {
          Interactions[i,1] <- substring(Source_reac[ReacIdx],2)
          Interactions[i,2] <- "-|"
        } else {
          Interactions[i,1] <- Source_reac[ReacIdx]
          Interactions[i,2] <- "->"
        }
        Interactions[i,3] <- Target_reac[ReacIdx]
        Interactions[i,4] <- paste('k',toString(counter+j),sep="")
        i <- i+1
      }
    }
  }

  # Add explicit inhibitory inteaction as an interaction
  if (!is.null(Inhibitor_names)) {
    for (counter in 1:length(Inhibitor_names)) {
      Interactions[i,1] <- Inhibitor_names_with_i[counter]
      Interactions[i,2] <- "-|"
      Interactions[i,3] <- Inhibitor_names[counter]
      Interactions[i,4] <- paste('ki',toString(counter),sep="")
      i <- i+1
    }
  }

  # Adding additional OR gate to the Interactions
  OR_idx <- which(Interactions[,5]=="O")
  if (length(OR_idx)>0) {
    OR_params <- unique(Interactions[OR_idx,4])
    OR_params_idx <- 1
    for (counter in 1:length(OR_params)) {
      Current_OR_idx <- which(Interactions[,4]==OR_params[counter])
      Current_OR_inputs <- Interactions[Current_OR_idx,1]
      Current_OR_output <- unique(Interactions[Current_OR_idx,3])
      Interactions[Current_OR_idx,3] <- paste(Current_OR_inputs[1],"_OR_",Current_OR_inputs[2],sep="")
      Interactions[Current_OR_idx,4] <- 1
      Reac1 <- c(Interactions[Current_OR_idx[1],1],"->", Current_OR_output,paste("kn",toString(OR_params_idx),sep=""),"N","D")
      Reac2 <- c(Interactions[Current_OR_idx[2],1],"->", Current_OR_output,paste("kn",toString(OR_params_idx+1),sep=""),"N","D")
      Reac3 <- c(paste(Interactions[Current_OR_idx[1],1],"_OR_",paste(Interactions[Current_OR_idx[2],1]),sep=""),"->", Current_OR_output,paste("kn",toString(OR_params_idx+2),sep=""),"N","D")
      Interactions <- rbind(Interactions,Reac1,Reac2,Reac3)
      state_names <- c(state_names,paste(Interactions[Current_OR_idx[1],1],"_OR_",paste(Interactions[Current_OR_idx[2],1]),sep=""))
      OR_params_idx <- OR_params_idx+3
    }
  }


  # Add manually assigned OR gate interaction
  if (!is.null(ORlist)) {
    ORlistNames <- matrix("NA",length(ORlist),3)
    for (counter in 1:(length(ORlist)/2)) {
      Splitted_SrcTar <- strsplit(ORlist[counter],split = "=")
      Splitted_ORchar <- strsplit(Splitted_SrcTar[[1]][1],split = "_OR_")
      ORlistNames[counter,1] <- Splitted_ORchar[[1]][1]
      ORlistNames[counter,2] <- Splitted_ORchar[[1]][2]
      ORlistNames[counter,3] <- Splitted_SrcTar[[1]][2]
    }
    for (counter in 1:nrow(ORlistNames)) {
      OR_1stIdx <- which(ORlistNames[counter,1]==Interactions[,1])
      OR_2ndIdx <- which(ORlistNames[counter,2]==Interactions[,1])
      OR_TarIdx <- which(ORlistNames[counter,3]==Interactions[,3])
      OR_Idx <- c(intersect(OR_1stIdx,OR_TarIdx),intersect(OR_2ndIdx,OR_TarIdx))
      if (Interactions[OR_Idx[1],2]==Interactions[OR_Idx[2],2]) { # make sure that the same sign is assigned for OR gate
        Interactions[OR_Idx,5] <- "O"
        Interactions[OR_Idx,4] <- "1"
      } else {
        stop("Please ensure that assigned OR interactions have the same type of interaction (activation/inhibition)")
      }
    }
  }

  # Check if there is any node with only inhibitory edge -> automatically add basal activity
  AllOutput <- unique(Interactions[,3])
  kb_idx <- 1
  for (counter in 1:length(AllOutput)) {
    OutAct <- which(AllOutput[counter]==Interactions[,3] & Interactions[,2]=="->")
    OutInb <- which(AllOutput[counter]==Interactions[,3] & Interactions[,2]=="-|")
    if (length(OutAct)==0 & length(OutInb)>0) { # if there is only inhibitory interaction to a single node -> add positive basal interaction
      # Interactions <- rbind(Interactions,c("basal","->",AllOutput[counter],paste0("kb",kb_idx),"N","D")) # Optimise basal
      Interactions <- rbind(Interactions,c("basal","->",AllOutput[counter],1,"N","D")) # NOT optimise basal
      if (kb_idx==1) {
        state_names <- c(state_names,"basal")
        Input_names <- c(Input_names,"basal")
        Input_index <- cbind(Input_index,rep(which(state_names=="basal"),nrow(Input_index)))
        Input_vector <- cbind(Input_vector,rep(1,nrow(Input_index)))
        colnames(Input_vector)[ncol(Input_vector)] <- "basal"
      }
      kb_idx <- kb_idx+1
    }
  }


  # If Forced, the weights of all single positive interactions are set to 1 (FALCON)
  if (Force==TRUE) {
    # Grab the interactions with all positive interactions
    Pos_IntAct <- which(Interactions[,2]=="->")
    Pos_Target <- Interactions[Pos_IntAct,3]
    Pos_Target_tab <- table(Pos_Target)
    Pos_Target_single <- names(which(Pos_Target_tab==1))
    for (counter in 1:length(Pos_IntAct)) {
      if (sum(Interactions[Pos_IntAct[counter],3]==Pos_Target_single)==1 || Interactions[Pos_IntAct[counter],5]=="A" || Interactions[Pos_IntAct[counter],5]=="O") {
        Interactions[Pos_IntAct[counter],4] <- "1"
      }
    }
  }

  Interactions_List <- list()
  for (counter in 1:nrow(Interactions)) {
    Interactions_List[[counter]] <- Interactions[counter,]
  }

  # ===== Step 3 : Build constraints from the list of interactions ==== #

  OptIn <- CNORprob_writeConstraint(Interactions,HLbound,state_names,HardConstraint)

  estim                   <<- list()
  estim$Interactions       <- Interactions
  estim$Interactions_List  <- Interactions_List
  estim$SSthresh           <- SSthresh
  estim$Input_vector       <- Input_vector
  estim$Input_index        <- Input_index
  estim$Output_vector      <- Output_vector
  estim$Output_index       <- Output_index
  estim$SD_vector          <- SD_vector
  estim$ma                 <- OptIn$ma
  estim$mi                 <- OptIn$mi
  estim$state_names        <- OptIn$state_names
  estim$param_vector       <- OptIn$param_vector
  estim$param_index        <- OptIn$param_index
  estim$FixBool            <- OptIn$FixBool
  estim$NrStates           <- length(OptIn$state_names)
  estim$NrParams           <- length(OptIn$param_vector)
  estim$PlotIterations     <- PlotIterations
  estim$rsolnp_options     <- rsolnp_options
  estim$A                  <- OptIn$A
  estim$b                  <- OptIn$b
  estim$Aeq                <- OptIn$Aeq
  estim$beq                <- OptIn$beq
  estim$LB                 <- OptIn$LB
  estim$UB                 <- OptIn$UB
  estim$L1Reg              <- L1Reg

  return(estim)

}

# ======= End of the script ======= #
