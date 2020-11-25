#' A function to run vatiance analysis
#'
#' @param contractid Character. The name of the contract or other value
#' @param mdata Data Frame. A filtered data frame of M-values
#' @param sampledata Data Frame. A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Character. Where to save the plots. (Default: working directory)
#' @return Outputs heatmaps and variance spreadsheets to your WD
#' @export
arrayvar <- function(contractid, mdata, sampledata, workdir=NULL) {
  for (package in c("RColorBrewer", "grDevices", "limma", "gplots")) {
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

  if(missing(sampledata)){
    message("error: sampledata is missing…")
    stop()
  }

  if(missing(mdata)){
    message("error: mdata is missing…")
    stop()
  }

  ## Set paths
  dir.create(paste0(workdir,"/secondary_analysis/Results/variance"), showWarnings = FALSE, recursive=TRUE)
  outprefix=paste0(workdir,"/secondary_analysis/Results/variance/", contractid)

  ## Set palette for MDS plots
  ngroups=ifelse(length(unique(sampledata$Sample_Group)) < 3, 3, length(unique(sampledata$Sample_Group)))
  pal <- RColorBrewer::brewer.pal(ngroups,"RdBu")
  col <- grDevices::colorRampPalette(pal)

  ## Set palette for heatmap
  pal2 <- RColorBrewer::brewer.pal(ngroups,"Set2")
  col.group <- RColorBrewer::brewer.pal(ngroups,"Set2")[as.factor(targets$Sample_Group)]

  ## Plot top variable genes
  vargenes <- apply(mdata, 1, var)
  vargenes500 <- names(sort(vargenes, decreasing=TRUE))[1:500]
  vargenes500_mdata <- mdata[vargenes500,]
  vargenes50 <- names(sort(vargenes, decreasing=TRUE))[1:50]
  vargenes50_mdata <- mdata[vargenes50,]

  ## Sources of variation analysis
  png(paste(outprefix,"_MDS.png", sep=""), res=150, height=700, width=1400)
    par(cex=0.4, mfrow=c(1,2))
    limma::plotMDS(mdata, top=1000, gene.selection = "common", col=pal2[factor(targets$Slide)], cex=0.6,
          main="Batch MDS - C1 v C2")
    legend("bottomleft",legend=levels(factor(targets$Slide)),text.col=pal2)
    limma::plotMDS(mdata, top=1000, gene.selection = "common", col=col.group, cex=0.6,
          main="Group MDS - C1 v C2")
    legend("bottomleft",legend=levels(factor(targets$Sample_Group)),text.col=pal2)
  dev.off()

  ## Examine higher dimensions to look at other sources of variation
  png(paste(outprefix, "_MDS_4dims.png", sep=""), res=150, height=500, width=1200)
    par(cex=0.5, mfrow=c(1,3))
    limma::plotMDS(mdata, top=1000, gene.selection = "common", col=col.group, dim.plot=c(1,3),
          cex=0.6, main="MDS - C1 v C3")
    legend("topright",legend=levels(factor(targets$Sample_Group)), text.col=pal2)
    limma::plotMDS(mdata, top=1000, gene.selection = "common", col=col.group, dim.plot=c(2,3),
          cex=0.6, main="MDS - C2 v C3")
    legend("topright",legend=levels(factor(targets$Sample_Group)), text.col=pal2)

    ## If there are enough samples then plot 3x4, else just plot top 3
    if (nrow(sampledata) > 4) {
      limma::plotMDS(mdata, top=1000, gene.selection = "common", col=col.group, dim.plot=c(3,4),
          cex=0.6, main="MDS - C3 v C4", ndim=nrow(targets)-1)
      legend("topright",legend=levels(factor(targets$Sample_Group)), text.col=pal2)
    }
  dev.off()

  png(paste(outprefix,"_PCA.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    plot(prcomp(mdata), main="Sources of variance", legend="bottomleft")
  dev.off()

  png(paste(outprefix,"_variance_heatmap_top500.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    gplots::heatmap.2(vargenes500_mdata, col=rev(col(50)), trace="none", main="Top 500 most variable genes across samples",
            ColSideColors=col.group,scale="row", labRow="")
    legend("topright",legend=levels(factor(targets$Sample_Group)), text.col=col.group, cex=0.5)
  dev.off()

  png(paste(outprefix,"_variance_heatmap_top50.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
    gplots::heatmap.2(vargenes50_mdata, col=rev(col(50)), trace="none", main="Top 50 most variable genes across samples",
            ColSideColors=col.group,scale="row")
    legend("topright",legend=levels(factor(targets$Sample_Group)), text.col=col.group, cex=0.5)
  dev.off()
}
