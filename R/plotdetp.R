#' Visualises the detection P value of array prior to filtering
#'
#' @param contractid Character. The name of the contract
#' @param pvals Data Frame. Detection P-values
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata()
#' @param workdir Character. Path to output location. (Default is working directory)
#' @return A barplot of the detection P values for the unfiltered array
#' @export
.plotDetP <- function(contractid, pvals, sampledata, workdir) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("The RColorBrewer package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(workdir)){
    message("error: workdir is missing......")
    exit()
  }

  if(missing(contractid)){
    message("error: contractid is missing......")
    exit()
  }

  if(missing(sampledata)){
    message("error: workdir is missing......")
    exit()
  }

  if(missing(pvals)){
    message("error: input 'pvals' not found.....")
    exit()
  }

  ## Take the mean
  plotdata = apply(pvals, 2 ,mean)

  ## Plot Detection P
  png(paste(workdir,"/secondary_analysis/Results/qc/", contractid,"_preQC_detectionP.png", sep=""), res=200, height=1000, width=1000)
    barplot(plotdata, las=2, cex.names=0.5, cex.axis=0.75,
        names.arg=sampledata$Sample_Name,
        col=pal[as.factor(sampledata$Sample_Group)],
        main="Mean detection p-values by sample",
        ylab="Mean detection p-value")
    abline(h=0.05, col="red")
    legend("bottomright", legend=levels(factor(sampledata$Sample_Group)), xpd=T, fill=pal, bg="white")
    try(dev.off(), silent = TRUE)
}
