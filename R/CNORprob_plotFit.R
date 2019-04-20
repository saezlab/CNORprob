#' @title CNORprob_plotFit
#'
#' @description Plot optimisation results against the experimental data
#'
#' @export

CNORprob_plotFit<-function(model, CNOlist, estim, res, show=TRUE, plotPDF=TRUE, tag=NULL, plotParams=list(cex=0.8, cmap_scale=1, yMin=0)) {

  if ((class(CNOlist)=="CNOlist")==FALSE){
    CNOlist = CellNOptR::CNOlist(CNOlist)
    # CNOlist = CNORuniv::CNOlist(CNOlist)
  }

  SimParams     <- res$BestFitParams
  SimIterations <- 3
  SimResAll     <- CNORprob_simulate(estim,SimParams,SimIterations)
  OutputIdx     <- estim$Output_index[1,]
  # simResults    <- list(t0=matrix(data=0,nrow=dim(SimResAll$MeanStateValueAll[,OutputIdx])[1],ncol=dim(SimResAll$MeanStateValueAll[,OutputIdx])[2]),t1=SimResAll$MeanStateValueAll[,OutputIdx])
  simResults    <- list(t0=CNOlist@signals[[1]],t1=SimResAll$MeanStateValueAll[,OutputIdx])

  if (show == TRUE){
    plotOptimResultsPan(
        simResults=simResults,
            CNOlist=CNOlist,
            formalism="ss1",
            tPt=CNOlist@timepoints[2],
            plotParams=plotParams)
  }

  if(plotPDF == TRUE){
    mkdirs("Results/Figures")
    if ( is.null(tag)){
           filename<-paste("Results/Figures/CNORprob_plotFit.pdf", sep="")
    } else {
        filename<-paste("Results/Figures/",tag, "CNORprob_plotFit.pdf", sep="_")
    }
    plotOptimResultsPan(
        simResults=simResults,
        CNOlist=CNOlist,
        formalism="ss1",
        tPt=CNOlist@timepoints[2],
        pdf=TRUE,
        pdfFileName=filename,
        plotParams=plotParams)
    }
  return(list(simResults=simResults))
}
