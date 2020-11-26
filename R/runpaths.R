#' Pathways Analysis
#'
#' A function to run pathways analysis from DM results
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrastsm Matrix. A matrix of contrasts
#' @param myrds List. An RDS file with a list of DM results
#' @param sampledata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Character. Where to save the plots. (Default: working directory)
#' @return CSV spreadsheets with Reactome pathways results and enrichment results. Dot plot and barplots of enrichment results.
#' @export
arraypaths <- function(contractid, sampledata, workdir=NULL, myrds, contrastm) {
  for (package in c("reactome.db", "biomaRt", "org.Hs.eg.db", "ReactomePA", "dplyr")) {
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
  dir.create(paste0(workdir,"/secondary_analysis/Results/pathways/"), showWarnings = FALSE)
  outprefix=paste0(workdir,"/secondary_analysis/Results/pathways/", contractid)

  ## Set inputs
  lrt_list <- readRDS(paste0(workdir, "/secondary_analysis/", "lrt_list.RDS"))
  comparisons <- colnames(contrastm)

  ## Load database
  print("Loading ensemble database......")
  ensembl=biomaRt::useMart("ensembl")
  ensembl=biomaRt::useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)
  gene_name_to_entrez = biomaRt::getBM(attributes=c("entrezgene_id", "external_gene_name"), mart=ensembl)

  i=1
  for (comparison in comparisons) {
    print(paste0("Running pathways for comparison: ", comparison, "......"))

    # Get names of samples used in comparison
    comparison.groups <- names(which(contrastm[,comparison] != 0))
    samples <- sampledata[sampledata$Sample_Group %in% comparison.groups, "Sample_Name"]

    # get comparison data for all probes
    #resall <- limma::topTable(fit2, adjust="BH", num=Inf, coef=comparison)
    resall <- lrt_list[[i]]

    # Annotate all probes with geneid
    .ann450k_tmp <- .ann450k[row.names(resall), ]
    genes <- strsplit(.ann450k_tmp$UCSC_RefGene_Name, ";", fixed=TRUE)
    genelabels <- sapply(genes, "[", 1)
    resall$genes <- genelabels

    # Get entrez vector in order of external gene names
    res_entrez_ids = gene_name_to_entrez[match(resall$genes, gene_name_to_entrez$external_gene_name), "entrezgene_id"]
    anno = data.frame(res_entrez_ids)
    rownames(anno) = rownames(resall)
    colnames(anno) = "Entrez ID"

    # Get entrezID to match to reactome pathway
    resall$`Entrez ID` <- anno$`Entrez ID`
    resall_pass <- resall[which(resall$`Entrez ID`!="NA"),]
    resall_anno <- resall_pass[ resall_pass$`Entrez ID` %in% AnnotationDbi::keys( reactome.db::reactome.db, "ENTREZID" ) & !is.na( resall_pass$adj.P.Val) , ]
    cat("nrows is: ")
    cat("ID's before FDR filter ", nrow(resall_pass),"\n")
    cat("ID's after filter ", nrow(resall_anno), "\n")

    ## Use all probes(not just unique genes)
    resall_anno$entrezid <- resall_anno$`Entrez ID`

    ## Annotate subset of DB (make unique keys)
    reactomeTable <- suppressWarnings(AnnotationDbi::select( reactome.db::reactome.db, keys=as.character(resall_anno$entrezid), keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID")))
    reactomeTable <- unique(reactomeTable)

    ## Convert Entrezt back to common gene name
    geneTable <- suppressWarnings(unique(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=reactomeTable$ENTREZID, keytype="ENTREZID", columns=c("SYMBOL","ENTREZID"))))
    reactomeDB_human<-merge(reactomeTable, geneTable, by="ENTREZID", all.x=T)[2:3]

    ## Make a vector of genes in each pathway
    reactomeDB_human_uniq <- unique(reactomeDB_human)
    res<-aggregate(SYMBOL~REACTOMEID, reactomeDB_human_uniq, FUN= paste, collapse='|')

    ## Final DF (matrix of reactome ID (rows) and entrezid (cols), TRUE/FALSE)
    reactDB <- do.call( rbind, with(reactomeTable, tapply(
      ENTREZID, factor(REACTOMEID), function(x) resall_anno$entrezid %in% x ) ))

    colnames(reactDB) <- rownames(resall_anno)

    ##remove paths with <20 genes or +80 genes
    within <- function(x, lower, upper) (x>=lower & x<=upper)
    reactDB <- reactDB[ within(rowSums(reactDB), lower=20, upper=80), ]

    ## Run t-test
    print("Running Pathways Analysis......")
    testCategory <- function( reactomeID ) {
      inReact <- reactDB[ reactomeID, ]
      myError<-tryCatch(t.test( resall_anno$logFC[ inReact ] )$p.value, error=function(e) e)
      data.frame(
        reactomeID = reactomeID,
        numGenes    = sum( inReact ),
        avgLFC      = mean( resall_anno$logFC[inReact] ),
        sdLFC       = sd( resall_anno$logFC[inReact] ),
        zValue      = mean( resall_anno$logFC[inReact] ) /sd( resall_anno$logFC[inReact] ),
        strength    = sum( resall_anno$logFC[inReact] ) / sqrt(sum(inReact)),
        pvalue=if(!inherits(myError, "error")){t.test( resall_anno$logFC[ inReact ] )$p.value} else {"NA"},
        #pvalue      = t.test( res2$logFC[ inReact ] )$p.value,
        reactomeName= if(is.character(reactome.db::reactomePATHID2NAME[[reactomeID]])){reactome.db::reactomePATHID2NAME[[reactomeID]]} else {print("NA")},
        #reactomeName = reactomePATHID2NAME[[reactomeID]][1],
        stringsAsFactors = FALSE, check.rows=T )
    }

    ## Build results
    tmp<-matrix(data = NA, ncol=8, nrow=nrow(reactDB)*2)
    k=1
    for (i in rownames(reactDB)) {
      tmp[k,] <- as.character(testCategory(i)[1,])
      k=k+1
    }

    results<-data.frame(tmp[1:nrow(reactDB), ])
    colnames(results) <- c("reactomeID","numGenes","avgLFC","sdLFC","zValue","strength","pvalue","reactomeName")
    colnames(res)<-c("reactomeID","symbol")
    results$padj <- p.adjust( as.numeric(as.character(results$pvalue)), "BH" )
    results<-dplyr::inner_join(results, res, by="reactomeID", all=T)
    pathway_results<-results[order(results$pvalue,decreasing = T), ]
    pathway_results_p05<-pathway_results[which(as.numeric(as.character(pathway_results$pvalue))<0.05), ]
    pathway_results_LFC<-head(results[order(abs(as.numeric(as.character(results$avgLFC))), decreasing=T), ], 50)

    write.table(pathway_results, paste(outprefix,"_", comparison, "_Reactome", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
    write.table(pathway_results_p05, paste(outprefix, "_", comparison, "_Reactome_Sig05", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
    write.table(pathway_results_LFC, paste(outprefix, "_", comparison, "_Reactome_top50_avgLFC", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")

    ## Run enrichment tests
    print("Running enrichment analysis......")
    resall_sig <- resall_anno[which(resall_anno$P.Value<0.05), ]
    n = ifelse(nrow(resall_sig)>500, 500, nrow(resall_sig))
    ngenes = length(unique(resall_sig$entrezid))

    if(ngenes>20) {
      print(paste0("Enrichment for top ", n, " probes"))
      enr <- ReactomePA::enrichPathway(gene=unique(resall_sig$entrezid), organism="human", pvalueCutoff=0.05, readable=T)

      png(paste0(outprefix, "_", comparison, "_enrichment_barplot_2000var.png"), res=125, width=1000, height=500)
        ReactomePA::barplot(enr, showCategory=10, title = paste0("barplot ", comparison))
      dev.off()

      png(paste0(outprefix, "_", comparison, "_enrichment_dotplot_2000var.png"), res=125, width=1000, height=500)
        ReactomePA::dotplot(enr, showCategory=10, title = paste0("dotplot ", comparison))
      dev.off()
    } else {
      print(paste0("Not enough significant genes for enrichemnt"))
    }
  }
}
