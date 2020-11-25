#' Plot QC results
#'
#' A function to run lumi QC
#'
#' @param mdata Data Frame. A filtered data frame of M-values
#' @param pdata Data Frame. A filtered data frame of detection p-values.
#' @param contractid Character. The name of the contract or other value
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata().
#' @param qctype Character indicating whether preQC or postQC is undertaken
#' @param workdir Character. Where to save the plots. (Default: working directory)
#' @return Default lumi QC plots
#' @export
.plotLumiQC <- function(contractid, mdata, sampledata, pdata, qctype, workdir=NULL) {
  for (package in c("lumi", "arrayQualityMetrics")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  if(missing(contractid)){
    message("error: contractid is missing…")
    stop()
  }

  if(missing(workdir)){
    workdir=getwd()
  }

  if(missing(mdata)){
    message("error: mdata are missing…")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing…")
    stop()
  }

  if(missing(pdata)){
    message("error: pdata is missing.....")
    stop()
  }

  if(missing(qctype)){
    message("error: qctype is missing.....")
    stop()
  }

  type <- ifelse(grep("pre", qctype, ignore.case=TRUE),
         grep("pre", qctype, value=TRUE, ignore.case=TRUE),
         grep("post", qctype, value=TRUE, ignore.case=TRUE))

  if (length(type) == 0) {
    message("error: type is invalid..... Input was: ", qctype)
    stop()
  }

  if(grep("pre", type, ignore.case=TRUE)){
    print("Running PreQC......")

    ## Set paths
    print("check3")
    dir.create(paste0(workdir,"/secondary_analysis/Results/qc"), showWarnings = FALSE, recursive=TRUE)
    outprefix=paste0(workdir,"/secondary_analysis/Results/qc/", contractid)

    ## Get quality stats, input to lumi
    colnames(mdata) <- sampledata$ID
    lumibatch <<- new('LumiBatch', exprs=mdata, detection=pdata, se.exprs=mdata)
    pData(lumibatch) <- sampledata
    lumistats <- arrayQualityMetrics::prepdata(expressionset = lumibatch, do.logtransform = FALSE, intgroup="Sample_Group")
    summary(lumibatch, 'QC')

  } else if(grep("post", type, ignore.case=TRUE)){
      print("Running PostQC......")

      ## Set paths
      outprefix=paste0(workdir,"/secondary_analysis/Results/qc/", contractid)

      ## Get quality stats, input to lumi
      lumibatchnorm <- new('LumiBatch', detection=pdata, exprs=mdata, se.exprs=mdata)
      pData(lumibatchnorm) <- sampledata
      lumistatsnorm <- arrayQualityMetrics::prepdata(expressionset = lumibatchnorm, do.logtransform = FALSE, intgroup="Sample_Group")
      summary(lumibatchnorm, 'QC')

      print("Plotting Gender......")
      png(paste(outprefix,"/../gender_check.png" ,sep=""), res=200, height=1000, width=1000)
        plotSex(rawquant)
      dev.off()

      png(paste(outprefix,"_", type, "_outliers.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
        detectOutlier(lumibatchnorm, ifPlot=TRUE)
      dev.off()

  } else {
      print("Error: qctype not value")
  }

  ## Plot heatmap
  print("OK1")
  png(paste(outprefix,"_", type, "_heatmap.png", sep=""), res=200, height=1000, width=1000)
   invisible(capture.output(arrayQualityMetrics::aqm.heatmap(lumistats)))
  dev.off()

  print("OK2")
  ## Generate QC plots
  png(paste(outprefix,"_", type, "_samplerelation.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    plot(lumibatch, what='sampleRelation')
  dev.off()

  print("OK3")
  png(paste(outprefix,"_", type, "_outliers.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    plot(lumibatch, what='outlier')
  dev.off()

  png(paste(outprefix,"_", type, "_samplerelation_MDS.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    plot(lumibatch, what='sampleRelation', method='mds')
  dev.off()

  png(paste(outprefix,"_", type, "_density_mvals.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    density(lumibatch, xlab="M-value")
  dev.off()

  png(paste(outprefix,"_", type, "_boxplot_mvals.png", sep=""), res=200, height=1000, width=1500); par(cex=0.6)
    boxplot(lumibatch)
  dev.off()

  png(paste(outprefix,"_", type, "_MDS_mvals.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    plotMDS(lumibatch)
  dev.off()
}

