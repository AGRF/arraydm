#' A function to run differential methylation
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrasts A named list with contrasts to apply
#' @param mdata A cleaned, quality controlled data frame of M-values
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Where to save the results (Default: working directory)
#' @return A data frame of modelling results
#' @export
arraymodel <- function(contractid, mdata, sampledata, contrastl, workdir=NULL) {
  for (package in c("limma")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  if(missing(contrastl)){
    message("error: contrastl are missing…")
    stop()
  }

  ## Set Paths
  dir.create(paste0(workdir,"/secondary_analysis/Results/comparisons"), showWarnings = FALSE, recursive=TRUE)
  outprefix=paste0(workdir,"/secondary_analysis/Results/comparisons/", contractid)

  ## Set up design and model
  print("Setting Contrasts......")
  design <- model.matrix(~0 + sampledata$Sample_Group)
  rownames(design)<-sampledata$Sample_Name
  colnames(design) <- gsub('sampledata\\$Sample_Group','', colnames(design))

  cont.matrix <<- eval(parse(text=paste0("limma::makeContrasts(",paste(names(contrastl),"=",contrastl, collapse=" ,"), ",levels=design)")))
  comparisons <- colnames(cont.matrix)
  print("Comparisons to run: ")
  print(paste0("    ",comparisons))
  print("Design matrix: ")
  print(paste0("    ", cont.matrix))
  write.table(cont.matrix, paste0(workdir, "/secondary_analysis/cont_matrix.txt"), quote=F)

  ## Linear modelling
  print("Running Differential Methylation Tests......")
  # or eset <- ExpressionSet(assayData = mdata)
  fit <- limma::lmFit(mdata, design)
  print("Results lmFit......")
  print(summary(limma::decideTests(fit)))

  ## Batch (SSM)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)

  ##Apply empirical Bayes smoothing to the standard errors.
  fit2 <- limma::eBayes(fit2)

  ## Summarise ebayes
  print("Results eBayes......")
  print(summary(limma::decideTests(fit2)))
  results <- limma::decideTests(fit2)

  ## Plot
  png(paste0(workdir, "/secondary_analysis/venn.png"), res=100, height=1000, width=1000)
    limma::vennDiagram(results)
  dev.off()

  ## Pairwise (SEM)
  lrt_list <- vector(mode="list", length=length(comparisons))
  names(lrt_list) = comparisons
  for (comparison in comparisons){
    fit2 <- limma::contrasts.fit(fit, cont.matrix[,comparison])
    fit2 <- limma::eBayes(fit2)

    # Save the table of results:
    print("OK1")
    lrt_list[[comparison]] <- limma::topTable(fit2, adjust="BH", number=Inf)
    print("OK2")
    res50 <- limma::topTable(fit2, adjust="BH", number=50, coef=comparison)
    print("OK3")
    resall <- limma::topTable(fit2, adjust="BH", num=Inf, coef=comparison)

    # write.table(res50, paste(outprefix,"_",comparison,"_LR_raw_top50.csv", sep=""), quote=F, sep="\t", col.names=T, row.names=T)
    write.table(resall, paste(outprefix,"_",comparison,"_LR_raw.csv", sep=""), quote=F, sep="\t", col.names=T, row.names=T)
  }

  ## For non-significant contrasts
#  for (comparison in comparisons){
#    resall <- topTable(fit2, adjust="BH", num=Inf, coef=comparison)
#    write.table(format(resall, digits=4), paste(outprefix,"_",comparison,"_LR_raw.csv", sep=""), quote=F, sep="\t", col.names=T, row.names=T)
#  }
}
