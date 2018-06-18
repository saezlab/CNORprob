#' @title CNORprob_BS
#'
#' @description Perform Bootstrapping with random
#'
#' @export

CNORprob_BS = function(model,CNOlist,estim,res,BS_Type,BS_Round) {

  estim_orig <<- estim
  model_orig <- model # Keep original CNOR model

  bestx <- res$BestFitParams
  num_params <- length(estim$param_vector)
  bestcost <- res$OptResults[1,1]
  L1Reg <- estim$L1Reg
  SSthresh <- estim$SSthresh

  cost_BS <- rep(0,1,BS_Round)
  param_BS <- matrix(NA,BS_Round,num_params)

  if (BS_Type==1) {
    ReorderedBS <- NULL
    for (counter in 1:BS_Round) {
      # set.seed(counter) # Optional
      ReorderedBS <- rbind(ReorderedBS, sample(x = 1:length(estim$Output_vector),size = length(estim$Output_vector), replace = T))
    }
    estim$ReorderedBS <- ReorderedBS

  } else if (BS_Type==2) {
    AddSD=0.05
    rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
    if (sum(estim$SD_vector,na.rm = T)==0) {
      estim$SD_vector <- matrix(rep(AddSD,length(estim$SD_vector)),nrow(estim$SD_vector),ncol(estim$SD_vector))
    }

    NewSampledData <- array(NA,c(dim(estim$Output_vector)[1],dim(estim$Output_vector)[2],BS_Round))
    for (counter in 1:length(estim$Output_vector)) {
      for (counter2 in 1:ncol(estim$Output_vector)) {
        for (counter3 in 1:nrow(estim$Output_vector)) {
          if (!is.na(estim$Output_vector[counter3,counter2])) {
            SampledData <- rnorm2(BS_Round,estim$Output_vector[counter3,counter2],estim$SD_vector[counter3,counter2])
          } else {
            SampledData <- rep(NA,BS_Round)
          }
          for (counter4 in 1:length(SampledData)) {
            NewSampledData[,,counter4][counter3,counter2] <- SampledData[counter4]
          }
        }
      }
    }

  }


  # for (counter in  1:length(p_KD)) {
  for (counter in  1:BS_Round) {

    print("=================================")
    print(paste("Bootstrapping Round:", toString(counter), "/", toString(BS_Round)))
    print("=================================")

    model_BS <- model_orig

    estim$maxtime <- estim_orig$maxtime
    estim$printCost <- estim_orig$printCost

    res <- NULL

    if (BS_Type==1) {
        estim$BS_Counter <- counter
    } else if (BS_Type==2) {
         estim$Output_vector <- NewSampledData[,,counter]
    }

    tryCatch(
      res <- CNORprob_optimise(estim,optRound=1)
    ,
    error=function(e) print("Solver failed... continue the next optimisation round"))

    if (!is.null(res)) {

      cost_BS[counter] <- res$OptResults[1,1] # fetch the best cost from optimizations
      param_BS[counter,] <- res$ParamResults[,2]

    } else { # if the optimiser failed
      cost_BS[counter] <- NA
      param_BS[counter,] <- rep(NA,num_params)
    }

  }

  # close(pb)

  # Plot figures
  estim <- estim_orig

  pdf("FitCost_Bootstrapping.pdf")
  boxplot(x = cost_BS, outpch = NA,main="Fitting Cost Bootstrapping")
  stripchart(x = cost_BS,
             vertical = TRUE, method = "jitter",
             pch = 21, col = "maroon", bg = "bisque",
             add = TRUE)
  dev.off()

  boxplot(x = cost_BS, outpch = NA,main="Fitting Cost Bootstrapping")
  stripchart(x = cost_BS,
             vertical = TRUE, method = "jitter",
             pch = 21, col = "maroon", bg = "bisque",
             add = TRUE)

  res$BestFitParams <- colMeans(param_BS)
  bString <- CNORprob_mapModel(model,estim,res)
  pdf("FitParam_Bootstrapping.pdf")
  plotModel(model=model,CNOlist = CNOlist,bString = bString)
  dev.off()
  plotModel(model=model,CNOlist = CNOlist,bString = bString)

  estim_Result$BS$Cost_BS      <- cost_BS
  estim_Result$BS$BS_param     <- bString
  estim_Result                <<- estim_Result
  estim                       <<- estim_orig
  estim                       <<- estim

  return(estim_Result)

}

# ======= End of the script ======= #
