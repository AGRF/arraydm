#' A function to collate and report differential methylation results
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrastm A matrix of contrasts
#' @param contrasts A named list with contrasts to apply
#' @param bdata A cleaned, quality controlled data frame of B-values
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Where to save the plots. (Default: working directory)
#' @return Annotated data frames of DM results, heatmaps and volcano plots
#' @export
arraysumm <- function(contractid, bdata, sampledata, workdir=NULL, contrastm) {
  for (package in c("RColorBrewer", "grDevices", "limma", "gplots")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  if(missing(contrastm)){
    message("error: contrastm is missing…")
    stop()
  }

  if(missing(contractid)){
    message("error: contractid is missing…")
    stop()
  }

  if(missing(workdir)){
    workdir=getwd()
  }

  if(missing(bdata)){
    message("error: bdata is missing…")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing…")
    stop()
  }

  ## Set outdir
  outprefix=paste0(workdir,"/secondary_analysis/Results/comparisons/", contractid)

  ## Annotate results
#  print("Annotating genes......")
#  genes <- strsplit(.ann450k$UCSC_RefGene_Name[1:50], ";", fixed=TRUE)
#  geneLabels <- sapply(genes, "[", 1)
#  genelabels <- matrix(geneLabels)
#  row.names(genelabels)<-rownames(.ann450k)

#  genelabels[which(is.na(genelabels)), 1] <- "NA"
  print(contrastm)
  comparisons <- colnames(contrastm)

  for (comparison in comparisons) {
    print(paste0("Summarising comparison: ", comparison, "......"))
    # Get names of samples used in comparison
    comparison.groups <- names(which(contrastm[,comparison] != 0))
    samples <- sampledata[sampledata$Sample_Group %in% comparison.groups, "Sample_Name"]

    # get comparison data
    resall <- limma::topTable(fit2, adjust="BH", num=Inf, coef=comparison)

    # annotate top 50 probes
    .ann450k_tmp50 <- .ann450k[row.names(resall)[1:50],]
    resannot_50 <- merge(resall[1:50,], .ann450k_tmp50, by="row.names", all.x=TRUE)
    topgenes <- .ann450k_tmp50$UCSC_RefGene_Name
    tmp <- sapply(strsplit(topgenes,split=";") , head, 1)
    genelab <- sapply(tmp, paste, collapse = ";")
    genelab[which(genelab=="")] <- "NA"

    # annotate all probes
    .ann450k_tmpall <- .ann450k[row.names(resall),]
    resannot_all <- merge(resall, .ann450k_tmpall, by="row.names", all.x=TRUE)

    # write files
    write.table(resannot_50, paste(outprefix,"_", comparison, "_DM_table_top50.csv", sep=""), sep=",", quote=F, col.names=T, row.names=F)
    write.table(resannot_all, paste(outprefix,"_", comparison, "_DM_table_all.csv", sep=""), sep=",", quote=F, col.names=T, row.names=F)

    ## Set palette for MDS plots
    ngroups=ifelse(length(unique(sampledata$Sample_Group)) < 3, 3, length(unique(sampledata$Sample_Group)))
    pal <- RColorBrewer::brewer.pal(11,"RdBu")
    col <- grDevices::colorRampPalette(pal)
    col.group <- RColorBrewer::brewer.pal(ngroups,"Set2")[as.factor(sampledata$Sample_Group)]

    ## Plot heatmaps
    png(paste(outprefix, "_", comparison, "_DM_heatmap_top50.png", sep=""), res=150, height=1000, width=1000); par(cex=0.4)
      gplots::heatmap.2(bdata[row.names(resall)[1:50], ], col=rev(col(50)), trace="none",
          main=paste("Top 50 most differentially methylated\ngenes between ", comparison, sep=""),
          ColSideColors=col.group,scale="row")
      legend("topright",legend=levels(factor(sampledata$Sample_Group)), text.col=unique(col.group), cex=0.6)
    dev.off()

    png(paste(outprefix, "_", comparison, "_DM_heatmap_top50_anno.png", sep=""), res=150, height=1000, width=1000); par(cex=0.4)
      gplots::heatmap.2(bdata[row.names(resall)[1:50], ], col=rev(col(50)), trace="none",
          main=paste("Top 50 most differentially methylated\ngenes between ", comparison, sep=""),
          ColSideColors=col.group,scale="row", labRow=genelab)
      legend("topright",legend=levels(factor(sampledata$Sample_Group)), text.col=unique(col.group), cex=0.6)
    dev.off()

    ## Pheatmap
    png(paste(outprefix, "_", comparison, "_DM_Pheatmap_top50.png", sep=""), res=150, height=1000, width=1000); par(cex=0.6)
      pheatmap::pheatmap(bdata[row.names(resall)[1:50],], color = viridis(100), trace="none",
         main=paste("Top 50 most differentially methylated\ngenes between ", comparison, sep=""), ColSideColors=col.group,
         scale="row", annotation_names_row=T, labels_row=genelab, fontsize_row=6)
    dev.off()

    ## Plot volcano
    .plotVolc(contractid, sampledata, workdir, contrastm)
  }
}
