#' A function to run gene ontology analysis from DM results
#'
#' @param contractid Character. The name of the contract or other value
#' @param myrds A list with DM results
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Where to save the plots. (Default: working directory)
#' @return CSV files with GO results
#' @export
dmgo <- function(contractid, sampledata, workdir=NULL, myrds) {
  for (package in c("GOstats", "lumiHumanAll.db", "DBI", "annotate")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  ## Run gene ontology for each comparison
  for (comparison in comparisons) {
    # Get names of samples used in comparison
    comparison.groups <- names(which(cont.matrix[,comparison] != 0))
    samples <- TARGETS[TARGETS$Sample_Group %in% comparison.groups, "Sample_Name"]

    # Get comparison data
    resall <- topTable(fit2, adjust="BH", num=Inf, coef=comparison)

    # Get significant genes only
    sigGene=resall[which(resall$P.Value<0.05),]

    # Add Gene IDs
    # Annotate all probes with geneid
    ann450k_tmp <- ann450k[row.names(sigGene),]
    allgenes <- strsplit(ann450k_tmp$UCSC_RefGene_Name, ";", fixed=TRUE)
    allgenelabels <- sapply(allgenes, "[", 1)

    # Add entrez IDs
    sigGene$genename<-allgenelabels
    res_entrez_ids = gene_name_to_entrez[match(sigGene$genename, gene_name_to_entrez$external_gene_name), "entrezgene_id"]
    anno = data.frame(res_entrez_ids)
    rownames(anno) = rownames(sigGene)
    colnames(anno) = "Entrez ID"
    sigGene$`Entrez ID` <- anno$`Entrez ID`

    # Run GO
    if (require(GOstats) & require(lumiHumanAll.db)) {
      ## Convert the probe Ids as  Entrez Ids and make them unique
      #sigLL <- unique(unlist(lookUp(sigGene,'lumiHumanAll.db','ENTREZID')))
      sigLL <-unique(unlist((sigGene$`Entrez ID`)))
      sigLL <- as.character(sigLL[!is.na(sigLL)])
      params <- new("GOHyperGParams",
                  geneIds= sigLL,
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

      ## Here only show the significant GO terms of BP (Molecular Function)
      ## For other categories, just follow the same procedure.
      sigGO.Term <- getGOTerm(sigGO.ID)[["BP"]]

      fdr<-data.frame(gGhyp.fdr[gGhyp.fdr < 0.01])
      fdr_pva<-cbind(fdr,data.frame(gGhyp.pv[rownames(fdr)]))

      gores<-cbind(sigGO.Term, fdr_pva)
      colnames(gores)<-c("sigGO.Term", "FDR", "PValue")
      gores$GO_ID <-rownames(gores)
      gores<-gores[ ,c(4,2,3,1)]
  }

  write.table(gores, paste(CONTRACT, "_", comparison, "_GO.csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
}
