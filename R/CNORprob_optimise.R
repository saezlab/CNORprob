#' @title CNORprob_optimise
#'
#' @description Optimisation interface between model descriptions 'estim' and the optimiser 'rsolnp'
#'
#' @export

CNORprob_optimise = function(estim,optRound=3,SaveOptResults=TRUE) {

  options(warn=-1)

  # Assign optimisation variables from estim
  MaxTime <- estim$maxtime

  # constraint functions
  A=estim$A; b=estim$b; Aeq=estim$Aeq; beq=estim$beq; LB=estim$LB; UB=estim$UB;

  # inequalities
  ineqn <- function(k) {A %*% k}
  ineqUB <- b

  # equalities
  eqn <- function( k ) {Aeq %*% k}
  eqB <- beq

  # Initialise variables to collect results
  resAll     <- list()
  ElapsedAll <- list()

  for (counter in 1:optRound) {

    k <- runif(length(estim$param_vector),min=0,max=1)
    while ((sum(k>UB) + sum(k<LB)) > 0) { # if initialised values are out of bound -> do it again
      k <- runif(length(estim$param_vector),min=0,max=1)
    }

    # Record the time (only for the optimisation process)
    ptm <- proc.time()
    res <- NULL

    if (length(A)==0 & length(Aeq)==0) {
      tryCatch({ # Use tryCatch to skip failed run(s)
        res <- withTimeout({
                  solnp( k,
                  CNORprob_optFunc,
                  LB=as.vector(LB),
                  UB=as.vector(UB),
                  control=estim$rsolnp_options)
                  }, timeout=MaxTime)
                },
                TimeoutException=function(ex) {
                cat("Timeout. Skipping.\n");
              },
              error=function(e) print("Solver failed... continue the next optimisation round"))
              # warning=function(w) print("Some warnings... continue the optimisation process"))
    } else if (length(A)==0 & length(Aeq)!=0) {
      tryCatch({ # Use tryCatch to skip failed run(s)
        res <- withTimeout({
                  solnp( k,
                  CNORprob_optFunc,
                  eqfun=eqn,
                  eqB=eqB,
                  LB=as.vector(LB),
                  UB=as.vector(UB),
                  control=estim$rsolnp_options)
                  }, timeout=MaxTime)
                },
                TimeoutException=function(ex) {
                cat("Timeout. Skipping.\n");
              },
              error=function(e) print("Solver failed... continue the next optimisation round"))
              # warning=function(w) print("Some warnings... continue the optimisation process"))
    } else if (length(A)!=0 & length(Aeq)==0) {
      tryCatch({ # Use tryCatch to skip failed run(s)
        res <- withTimeout({
                  solnp( k,
                  CNORprob_optFunc,
                  ineqfun=ineqn,
                  ineqLB=rep(0,length(ineqUB)),
                  ineqUB=ineqUB,
                  LB=as.vector(LB),
                  UB=as.vector(UB),
                  control=estim$rsolnp_options)
                }, timeout=MaxTime)
                },
                TimeoutException=function(ex) {
                cat("Timeout. Skipping.\n");
              },
              error=function(e) print("Solver failed... continue the next optimisation round"))
              # warning=function(w) print("Some warnings... continue the optimisation process"))
    } else {
      tryCatch({ # Use tryCatch to skip failed run(s)
        res <- withTimeout({
                  solnp( k,
                  CNORprob_optFunc,
                  ineqfun=ineqn,
                  ineqLB=rep(0,length(ineqUB)),
                  ineqUB=ineqUB,
                  eqfun=eqn,
                  eqB=eqB,
                  LB=as.vector(LB),
                  UB=as.vector(UB),
                  control=estim$rsolnp_options)
                }, timeout=MaxTime)
                },
                TimeoutException=function(ex) {
                cat("Timeout. Skipping.\n");
              },
              error=function(e) print("Solver failed... continue the next optimisation round"))
              # warning=function(w) print("Some warnings... continue the optimisation process"))
    }

    if (exists("res")) {
      print( res$pars )
      resAll[[counter]]     <- res
      ElapsedAll[[counter]] <- proc.time() - ptm
    } else {
      res <- NULL
    }
  }


  if (exists("res")) {

    # Extract fitting costs and optimised paramters
    FittingCostAll <- rep(NaN,length(resAll))
    OptParamsAll   <- matrix(NaN,length(resAll),length(estim$param_vector))
    for (counter in 1:length(resAll)) {
      if (length(resAll[[counter]])>0) {
        FittingCostAll[counter] <- tail(resAll[[counter]]$values,n=1)
        OptParamsAll[counter,]   <- resAll[[counter]]$pars
      }
    }

    # Remove empty results from failed optimisation runs [i.e. remove NaN if any]
    FittingCostAll <- FittingCostAll[!is.na(FittingCostAll)]
    OptParamsAll <- OptParamsAll[!is.na(OptParamsAll)]
    # dim(OptParamsAll) <- c(length(OptParamsAll)/length(OptIn$param_vector),length(OptIn$param_vector))
    dim(OptParamsAll) <- c(length(OptParamsAll)/length(estim$param_vector),length(estim$param_vector))

    # Find the bestrun (lowerest fitting cost) and extract the corresponding parameter values
    BestFitIndex  <- which(FittingCostAll==min(FittingCostAll))
    BestFitParams <- OptParamsAll[BestFitIndex[1],] # Note, if multiple indices, accept the first one

    # Report optimisation results
    OptResults <- as.data.frame(matrix(nrow = 2,ncol = 3))
    colnames(OptResults) <- c("best","mean","S.D.")
    rownames(OptResults) <- c("FitCost","Time(s)")
    OptResults[1,1] <- min(FittingCostAll)
    OptResults[1,2] <- mean(FittingCostAll)
    OptResults[1,3] <- sd(FittingCostAll)

    ElapsedTime <- NULL
    for (ElapsedCount in 1:length(ElapsedAll)) {
      ElapsedTime <- c(ElapsedTime,ElapsedAll[[ElapsedCount]][3])
    }

    OptResults[2,1] <- min(ElapsedTime)
    OptResults[2,2] <- mean(ElapsedTime)
    OptResults[2,3] <- sd(ElapsedTime)

    print("Summarised fitting cost and time:")
    print(OptResults)

    # Report fitted parameters
    ParamResults <- as.data.frame(matrix(data = NA,nrow = length(estim$param_vector),ncol = 4))
    colnames(ParamResults) <- c("parameters","best","mean","S.D.")
    ParamResults[,1] <- as.character(estim$param_vector)
    ParamResults[,2] <- round(BestFitParams,digits=4)
    ParamResults[,3] <- round(colMeans(OptParamsAll),digits=4)
    ParamResults[,4] <- round(apply(OptParamsAll,2,sd),digits=4)

    print("Summarised identified parameters:")
    print(ParamResults)

    # Save results
    if (SaveOptResults) {
      mkdirs("Results/Optimisation")
      write.table(OptResults,file="Results/Optimisation/Optimisation_Results_CNORprob.txt",quote=FALSE,row.names=TRUE,col.names = TRUE)
      write.table(ParamResults,file="Results/Optimisation/Optimised_Parameters_CNORprob.txt",quote=FALSE,row.names=FALSE,col.names = TRUE)
      save(resAll,file='Results/Optimisation/OptResults_CNORprob.RData')
      save(ElapsedAll,file='Results/Optimisation/OptTime_CNORprob.Rdata')
    }

    res <- list()
    res$OptResults <- OptResults
    res$ParamResults <- ParamResults
    res$resAll <- resAll
    res$ElapsedAll <- ElapsedAll
    res$BestFitParams <- BestFitParams

  } else {
    res <- NULL
  }

  estim_Result$optimisation <<- res

  options(warn=0)

  return(res)

}
