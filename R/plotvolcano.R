#' A function to generate volcano plots from DM results
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrastm Matrix. A matrix of contrasts
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Charater. Where to save the plots. (Default: working directory)
#' @return One volcano plot per comparison
#' @export
.plotVolc <- function(contractid, sampledata, workdir=NULL, contrastm) {
  for (package in c("RColorBrewer", "grDevices", "limma", "gplots")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  if(missing(workdir)){
    workdir=getwd()
  }

  if(missing(contractid)){
    message("error: contractid is missing......")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing......")
    stop()
  }

  if(missing(contrastm)){
    message("error: contrastm is missing…")
    stop()
  }

  ## Set outdir
  outprefix=paste0(workdir,"/secondary_analysis/Results/comparisons/", contractid)

  ## Set inputs
  lrt_list <- readRDS(paste0(workdir, "/secondary_analysis/", "lrt_list.RDS"))
  comparisons <- colnames(contrastm)

  i=1
  ## Make a volcano plot for each comparison
  for (comparison in comparisons) {
    # Get names of samples used in comparison
    comparison.groups <- names(which(contrastm[,comparison] != 0))
    samples <- sampledata[sampledata$Sample_Group %in% comparison.groups, "Sample_Name"]

    # get comparison data
    #resall <- limma::topTable(fit2, adjust="BH", num=Inf, coef=comparison)
    resall <- lrt_list[[i]]

    ## palette is Set2 from colorbrewer
    green="#66c2a5"
    orange="#fc8d62"
    blue="#8da0cb"
    pink="#e78ac3"

    # annotate top probes "Green"
    resall_top <- subset(resall, adj.P.Val<.05 & abs(logFC)>1)
    .ann450k_tmp <- .ann450k[row.names(resall_top),]
    topgenes <- .ann450k_tmp$UCSC_RefGene_Name
    tmp <- sapply(strsplit(topgenes,split=";") , head, 1)
    genelab <- sapply(tmp, paste, collapse = ";")
    genelab[which(genelab=="")] <- "NA"

    # annotate all probes
    .ann450k_tmpall <- .ann450k[row.names(resall),]
    resannot_all <- merge(resall, .ann450k_tmpall, by="row.names", all.x=TRUE)

    png(paste(outprefix, "_", comparison, "_DM_volcano_plot.png", sep=""), res=200, height=1000, width=1000); par(cex=0.6)
      with(resall, plot(logFC, -log10(P.Value), pch=20, main=paste("Volcano plot\n",comparison, sep=""), xlim=c(-5,5)))

      # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
      with(subset(resall, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col=blue))
      with(subset(resall, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col=orange))
      with(subset(resall, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col=green))

      # Label points with the textxy function from the calibrate plot
      # cannot do if too many DM genes
      #	with(subset(resall, adj.P.Val<.05 & abs(logFC)>1), textxy(logFC, -log10(P.Value), labs=genelab, cex=.8))
    dev.off()

    i=i+1
  }
}
