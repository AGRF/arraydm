#' Gene Ontology
#'
#' A function to run gene ontology analysis from DM results
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrastsm Matrix. A matrix of contrasts
#' @param myrds List. An RDS file with lists of DM results
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Character. Where to save the plots. (Default: working directory)
#' @return CSV spreadsheets with Gene Ontology results.
#' @export
dmgo <- function(contractid, sampledata, workdir=NULL, myrds, contrastm) {
  for (package in c("GOstats", "biomaRt", "lumiHumanAll.db", "DBI", "annotate")) {
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

  if(missing(myrds)){
    message("error: myrds value is missing…")
    stop()
  }

  if(missing(contrastm)){
    message("error: contrastm is missing…")
    stop()
  }

  ## Set paths
  dir.create(paste0(workdir,"/secondary_analysis/Results/pathways"), showWarnings = FALSE)
  outprefix=paste0(workdir,"/secondary_analysis/Results/pathways/", contractid)

  ## Set inputs
  lrt_list <- readRDS(paste0(workdir, "/secondary_analysis/", "lrt_list.RDS"))
  comparisons <- colnames(contrastm)

  ## Load database
  print("Loading ensemble database......")
  ensembl=biomaRt::useMart("ensembl")
  ensembl=biomaRt::useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)
  gene_name_to_entrez = biomaRt::getBM(attributes=c("entrezgene_id", "external_gene_name"), mart=ensembl)

  ## Run gene ontology for each comparison
  i=1
  for (comparison in comparisons) {
    # Get names of samples used in comparison
    comparison.groups <- names(which(cont.matrix[,comparison] != 0))
    samples <- sampledata[sampledata$Sample_Group %in% comparison.groups, "Sample_Name"]

    # Get comparison data
    # resall <- topTable(fit2, adjust="BH", num=Inf, coef=comparison)
    resall <- lrt_list[[i]]

    # Get significant genes only
    resall_sig=resall[which(resall$P.Value<0.05),]

    # Add Gene IDs
    # Annotate all probes with geneid
    .ann450k_tmp <- .ann450k[row.names(resall_sig),]
    genes <- strsplit(.ann450k_tmp$UCSC_RefGene_Name, ";", fixed=TRUE)
    genelabels <- sapply(genes, "[", 1)

    # Add entrez IDs
    resall_sig$genes <- genelabels
    res_entrez_ids <- gene_name_to_entrez[match(resall_sig$genes, gene_name_to_entrez$external_gene_name), "entrezgene_id"]
    anno = data.frame(res_entrez_ids)
    rownames(anno) = rownames(resall_sig)
    colnames(anno) = "Entrez ID"
    resall_sig$`Entrez ID` <- anno$`Entrez ID`

    # Run GO
    if (require(GOstats) & require(lumiHumanAll.db)) {
      ## Convert the probe Ids as  Entrez Ids and make them unique
      geneids <- unique(unlist((resall_sig$`Entrez ID`)))
      geneids <- as.character(geneids[!is.na(geneids)])
      params <- new("GOHyperGParams",
                  geneIds= geneids,
                  annotation="lumiHumanAll.db",
                  ontology="BP",
                  pvalueCutoff= 0.01,
                  conditional=FALSE,
                  testDirection="over")
      hgOver <- hyperGTest(params)

      ## Get the p-values of the test
      gGhyp.pv <- pvalues(hgOver)

      ## Adjust p-values for multiple test (FDR)
      gGhyp.fdr <- p.adjust(gGhyp.pv,'fdr')

      ## select the Go terms with adjusted p-value less than 0.01
      sigGO.ID <- names(gGhyp.fdr[gGhyp.fdr < 0.01])

      ## Here only show the significant GO terms of BP, MF (Molecular Function)
      ## For other categories, just follow the same procedure.
      sigGO.Term <- annotate::getGOTerm(sigGO.ID)[["BP"]]

      ## Generate outputs
      fdr<-data.frame(gGhyp.fdr[gGhyp.fdr < 0.01])
      fdr_pva<-cbind(fdr,data.frame(gGhyp.pv[rownames(fdr)]))
      gores<-cbind(sigGO.Term, fdr_pva)
      colnames(gores)<-c("sigGO.Term", "FDR", "PValue")
      gores$GO_ID <-rownames(gores)
      gores<-gores[ ,c(4,2,3,1)]
   }

  write.table(gores, paste(outprefix, "_", comparison, "_GO.csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
  i=i+1
  }
}
