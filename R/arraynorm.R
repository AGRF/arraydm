#' Normalises methylation array data
#'
#' @param arraydata An RGChannelSet generated manually, or in arraydm::readdata()
#' @param contractid The name of the contract
#' @param mdataraw A filtered data frame of M-values
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata()
#' @param workdir Where to save the plots. (Default: working directory)
#' @return PreQC plots for detection p-value, quality filtered probe set
#' @export
arraynorm <- function(contractid, arraydata, mdataraw, sampledata, workdir=NULL) {
  if (!requireNamespace("arrayQualityMetrics", quietly = TRUE)) {
    stop("RColorBrewer is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(workdir)){
    workdir=getwd()
  }

  if(missing(contractid)){
    message("error: contractid is missing......")
    stop()
  }

  outprefix=paste0(workdir, "/secondary_analysis/Results/", contractid)

  if(missing(mdataraw)){
    message("error: mdataraw is missing…")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing…")
    stop()
  }

  if(missing(arraydata)){
    message("error: arraydata is missing…")
    stop()
  }

  ## Normalise
  print("Normalising Arrays..........")
  rawswan = preprocessSWAN(arraydata)
  rawquant = preprocessQuantile(arraydata)

  ## Plot Normalised mvals
  print("Plotting Results..........")
  .plotNorm(contractid, rawdata=mdataraw, quantdata=rawquant, swandata=rawswan, sampledata=targets, workdir)

  ## Swan normalise raw data GenomicMethylSet
  print("Proceeding with SWAN normalisation..........")
  swannorm = preprocessSWAN(rgSet = arraydata, mSet = mdataraw, verbose=TRUE)
  swannorm <<- mapToGenome(swannorm)
}
