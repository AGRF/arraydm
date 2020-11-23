#' Visualises the detection P value of array prior to filtering
#'
#' @param contractid The name of the contract
#' @param pvals Data Frame with detection P-values
#' @return A barplot of the detection P values for the unfiltered array
#' @export
.plotDetP <- function(contractid, pvals, workdir) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("The RColorBrewer package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(contractid)){
    message("error: contractid is missing......")
    exit()
  }

  if(missing(pvals)){
    message("error: input 'pvals' not found.....")
    exit()
  }

  dir.create(paste0(workdir,"/secondary_analysis/Results/preqc"), showWarnings = FALSE, recursive=TRUE)
  png(paste(workdir,"/secondary_analysis/Results/preqc/", contractid,"_preQC_detectionP.png", sep=""), res=200, height=1000, width=1000)
    barplot(apply(pvals, 2 ,mean), las=2, cex.names=0.5, cex.axis=0.75,
        names.arg=targets$Sample_Name,
        col=pal[as.factor(targets$Sample_Group)],
        main="Mean detection p-values by sample", ylab="Mean detection p-value")
    abline(h=0.05, col="red")
    legend("bottomright", legend=levels(factor(targets$Sample_Group)), xpd=T, fill=pal, bg="white")
  dev.off()
}
