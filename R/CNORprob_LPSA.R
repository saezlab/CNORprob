#' @title CNORprob_LPSA
#'
#' @description Systematically run local parameter sensitivity analysis 'LPSA' for each optimised parameter to assess their identifiabilities in the range of 0 to 1
#'
#' @export

CNORprob_LPSA = function(model,CNOlist,estim,res,HLbound=0.5,LPSA_Increments=2,Force=FALSE) {

  optRound_LPSA=estim$optRound_analysis
  LPSA_Increments=LPSA_Increments
  if (model$reacID[1] == "EGF=PI3K") { # if CNOToy model, then use ExpandOR section
    Force <- FALSE
    HardConstraint <- FALSE
    ExpandOR <- TRUE
  } else {
    Force <- TRUE
    HardConstraint <- TRUE
    ExpandOR <- FALSE
  }

  estim_orig <<- estim
  A <- estim$A
  Aeq <- estim$Aeq

  bestx <- res$BestFitParams
  SSthresh <- estim$SSthresh

  p <- bestx
  p_SA <- matrix(NA,2*LPSA_Increments+1,length(p))
  p_SA[LPSA_Increments+1,] <- p

  for (counter in 1:length(p)) {
    Min <- estim$LB[counter]
    Max <- estim$UB[counter]

    if (length(A)>0) {
      if (sum(A[,counter])>0) {
        IdxC <- which(A[,counter]>0);
        Max  <- min(estim$b[IdxC],Max);
      }
    }

    for (counter2 in 1:(LPSA_Increments-1)) {
      p_pert_values_pos <- p[counter] + ((Max-p[counter])/LPSA_Increments*counter2)
      p_pert_values_neg <- p[counter]-((p[counter]-Min)/LPSA_Increments*counter2)
      p_SA[LPSA_Increments+1+counter2,counter] = p_pert_values_pos
      p_SA[LPSA_Increments+1-counter2,counter] = p_pert_values_neg
    }

    p_SA[dim(p_SA)[1],counter] <- Max
    p_SA[1,counter] <- Min
  }

  # Pick index and evaluate
  cost_SA   <- matrix(NaN,dim(p_SA)[1],dim(p_SA)[2])
  params_SA <- list()
  for (counter in 1:dim(p_SA)[2]) {
    params_SA[[counter]] <- matrix(NaN,dim(p_SA)[1],dim(p_SA)[2])
  }

  # Resimulation with perturbed parameter values
  estim                <<- estim_orig
  Param_original        <- estim$param_vector
  State_names_original  <- estim$state_names
  Interactions_original <- estim$Interactions

  MidPoint              <- ceiling(dim(p_SA)[1]/2)

  num_plots             <- length(Param_original)
  NLines                <- ceiling(sqrt(num_plots))
  NCols                 <- ceiling(num_plots/NLines)

  Ident_All             <- c()

  # pb <- tkProgressBar(title = "LPSA analysis", min = 0, max = length(Param_original), width=300)

  estim$Input_vector    -> Input_vector
  estim$Input_index     -> Input_index
  estim$Output_vector   -> Output_vector
  estim$Output_index    -> Output_index
  estim$SD_vector       -> SD_vector

  for (counter in 1:dim(p_SA)[2]) { # for each parameter

    # setTkProgressBar(pb, counter, label=paste( round(counter/length(Param_original)*100, 0), "% done"))
    # setTkProgressBar(pb, counter, label=paste("Running LPSA on parameter:", toString(counter), "/", toString(length(Param_original))))

    print("===================================")
    print(paste("Running LPSA on parameter:", toString(counter), "/", toString(length(Param_original))))
    print("===================================")

    for (counter2 in 1:dim(p_SA)[1]) {

      Interactions <- Interactions_original

      Interactions_params <- Interactions[,4]
      replace_idx  <- which(Interactions_params==Param_original[[counter]])
      for (counter3 in 1:length(replace_idx)) {
        Interactions[replace_idx[counter3],4] <- toString(p_SA[counter2, counter])
      }

      if (Force==TRUE) {
        # Grab the interactions with all positive interactions
        Pos_IntAct <- which(Interactions[,2]=="->")
        Pos_Target <- Interactions[Pos_IntAct,3]
        Pos_Target_tab <- table(Pos_Target)
        Pos_Target_single <- names(which(Pos_Target_tab==1))
        for (counter4 in 1:length(Pos_IntAct)) {
          if (sum(Interactions[Pos_IntAct[counter4],3]==Pos_Target_single)==1 || Interactions[Pos_IntAct[counter4],5]=="A" || Interactions[Pos_IntAct[counter4],5]=="O") {
            Interactions[Pos_IntAct[counter4],4] <- "1"
          }
        }
      }

      state_names <- State_names_original

      OptIn <- CNORprob_writeConstraint(Interactions,HLbound,state_names,HardConstraint,LPSA=replace_idx)

      # Store all necessary elements into "estim" global variable

      L1Reg <- estim_orig$L1Reg
      SSthresh <- estim_orig$SSthresh
      PlotIterations <- estim_orig$PlotIterations

      Interactions <- Interactions[-replace_idx,]

      if (!is.null(dim(Interactions))) {
        Interactions_List <- list()
        for (counter_List in 1:nrow(Interactions)) {
          Interactions_List[[counter_List]] <- Interactions[counter_List,]
        }
      }

      estim                <<- list()
      estim$Interactions_List <- Interactions_List
      estim$Interactions    <- Interactions
      estim$SSthresh        <- SSthresh
      estim$Input_vector    <- Input_vector
      estim$Input_index     <- Input_index
      estim$Output_vector   <- Output_vector
      estim$Output_index    <- Output_index
      estim$SD_vector       <- SD_vector
      estim$ma              <- OptIn$ma
      estim$mi              <- OptIn$mi
      estim$state_names     <- OptIn$state_names
      estim$param_vector    <- OptIn$param_vector
      estim$param_index     <- OptIn$param_index
      estim$FixBool         <- OptIn$FixBool
      estim$NrStates        <- length(OptIn$state_names)
      estim$NrParams        <- length(OptIn$param_vector)
      estim$PlotIterations  <- PlotIterations
      estim$rsolnp_options  <- estim_orig$rsolnp_options
      estim$A               <- OptIn$A
      estim$b               <- OptIn$b
      estim$Aeq             <- OptIn$Aeq
      estim$beq             <- OptIn$beq
      estim$LB              <- OptIn$LB
      estim$UB              <- OptIn$UB
      estim$L1Reg           <- L1Reg
      estim$maxtime         <- estim_orig$maxtime
      estim$printCost       <- estim_orig$printCost

      resAll     <- list()
      ElapsedAll <- list()
      cost_LPSA  <- NULL

      for (counter_Round in 1:optRound_LPSA) {
        ptm <- proc.time()

        tryCatch(
        res <- CNORprob_optimise(estim,optRound=1,SaveOptResults = FALSE)
        ,
        error=function(e) print("Solver failed... continue the next optimisation round"))

        if (exists("res")) {
          # print( res$ParamResults )
          cost_LPSA <- rbind(cost_LPSA,c(counter_Round,res$OptResults[1,1]))
          resAll[[counter_Round]]     <- res
          ElapsedAll[[counter_Round]] <- proc.time() - ptm
        } else {
          cost_LPSA <- rbind(cost_LPSA,c(counter_Round,NA))
          resAll[[counter_Round]]     <- NA
          ElapsedAll[[counter_Round]] <- NA
        }
      }

      if (length(resAll)>0) {
        # cost_SA[counter2,counter] <- tail(resAll[[1]]$values,n = 1)
        Idx_BestCostLPSA <- cost_LPSA[which(cost_LPSA[,2]==min(cost_LPSA[,2])),1]
        cost_SA[counter2,counter] <- resAll[[Idx_BestCostLPSA[1]]]$OptResults[1,1]
        included_params_idx <- NULL
        for (counter7 in 1:length(estim$param_vector)) {
          included_params_idx <- c(included_params_idx, which(Param_original==estim$param_vector[[counter7]]))
        }
        for (counter8 in 1:length(included_params_idx)) {
          params_SA[[included_params_idx[counter8]]][counter2,counter] <- resAll[[Idx_BestCostLPSA[1]]]$ParamResults[counter8,2]
        }
      }
    }
  }

  print("New fitting cost of all parameters")
  print(cost_SA)
  print("Range of explore parameter values")
  print(params_SA)
  estim_Result$LPSA$cost_SA          <- cost_SA
  estim_Result$LPSA$params_SA        <- params_SA
  estim_Result$LPSA$p_SA             <- p_SA
  estim_Result$LPSA$param_vector     <- estim$param_vector
  estim_Result                       <<- estim_Result
  estim                              <<- estim_orig

  # Start ploting
  TotalParams <- dim(p_SA)[2]
  NLines=ceiling(sqrt(TotalParams))
  NCols=ceiling(TotalParams/NLines)
  par(mfrow=c(NLines,NCols),oma=c(0,0,2,0))

  pdf("Results/Figures/CNORprob_LPSA.pdf")
  for (counter in 1:TotalParams) {
    plot(p_SA[,counter],cost_SA[,counter],xlab="Parameter range",ylab="MSE",main=estim$param_vector[counter],type="o",col="blue",pch=16,cex=1.5)
    # if (!is.null(estim$SD_vector)) { # if variance of the data was provided
    #   lines(p_SA[,counter],rep(CutOff,dim(p_SA)[1]),col="red")
    # }
    points(p_SA[(dim(p_SA)[1]+1)/2,counter],cost_SA[(dim(p_SA)[1]+1)/2,counter],col="blue",pch=1,cex=2)
  }
  mtext("Local Parameters Sensitivity Analysis (LPSA)",outer=T,cex=1)
  dev.off()

  par(mfrow=c(1,1)) # Return the original full plot configuration for the next plot(s)

  print("LPSA results are saved as a figure in the Results folder")

  return(estim_Result)

}

# ======= End of the script ======= #
