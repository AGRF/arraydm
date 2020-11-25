#' Array Normalisation
#'
#' Plots normalised data. Best applied after arraydm::arraypreqc() where most of the inputs are generated.
#'
#' @param contract The name of the contract or similar ID
#' @param rawdata Data Frame. Raw M-values
#' @param quantdata Quantile normalised array.
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata()
#' @param swandata SWAN normalised array.
#' @param workdir Character. Where to save the plots. (Default: working directory)
#' @return Default SWAN Normalised array and data distribution
#' @export
.plotNorm <- function(contractid, rawdata, quantdata, swandata, sampledata, workdir){
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("RColorBrewer is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("minfi is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(contractid)){
    message("error: contractid is missing…")
    stop()
  }

  if(missing(workdir)){
    workdir=getwd()
  }

  if(missing(rawdata)){
    message("error: rawdata is missing…")
    stop()
  }

  if(missing(quantdata)){
    message("error: quantdata is missing…")
    stop()
  }

  if(missing(swandata)){
    message("error: swandata is missing…")
    stop()
  }

  ## Set Paths
  outprefix=paste0(workdir,"/secondary_analysis/Results/qc/", contractid)
  dir.create(paste0(workdir,"/secondary_analysis/Results/qc/"), showWarnings = FALSE, recursive=T)

  ## Set plotting palette
  ngroups=ifelse(length(unique(sampledata$Sample_Group)) < 3, 3, length(unique(sampledata$Sample_Group)))
  pal <<- RColorBrewer::brewer.pal(ngroups,"Set2")    # Bad practice but im doing it anyway

  ## Plot normalisations
  png(paste(outprefix,"_postQC_normalisation.png", sep=""), res=150, height=700, width=1000)
    par(mfrow=c(1,3))
    par(cex=0.6)
    minfi::densityPlot(rawdata, sampGroups=sampledata$Sample_Group,main="Raw", legend=FALSE)
    legend("top", legend = levels(factor(sampledata$Sample_Group)), text.col=RColorBrewer::brewer.pal(ngroups,"Set2"))

    minfi::densityPlot(getBeta(quantdata), sampGroups=sampledata$Sample_Group, main="Quantile Normalisation", legend=FALSE)
    legend("top", legend = levels(factor(sampledata$Sample_Group)), text.col=RColorBrewer::brewer.pal(ngroups,"Set2"))

    minfi::densityPlot(getBeta(swandata), sampGroups=sampledata$Sample_Group, main="SWAN Normalisation", legend=FALSE)
    legend("top", legend = levels(factor(sampledata$Sample_Group)), text.col=RColorBrewer::brewer.pal(ngroups,"Set2"))
  dev.off()
}
