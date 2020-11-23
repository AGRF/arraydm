#' A function to run pathways analysis from DM results
#'
#' @param contractid Character. The name of the contract or other value
#' @param myrds An RDS of DM results
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Where to save the plots. (Default: working directory)
#' @return CSV spreadsheets with Reactome pathways results
#' @export
dmpaths <- function(contractid, sampledata, workdir=NULL, myrds) {
  for (package in c("reactome.db", "biomaRt", "org.Hs.eg.db")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop(paste0("The ", package, " package is needed for this function to work. Please install it."),
           call. = FALSE)
    }
  }

  ## Load database
  print("Loading biomaRt......")
  ensembl=useMart("ensembl")
  ensembl=useDataset(dataset="hsapiens_gene_ensembl", mart=ensembl)
  gene_name_to_entrez = getBM(attributes=c("entrezgene_id", "external_gene_name"), mart=ensembl)


  ## load pre-processed data
  # rawswan_filt2.rds <- readRDS("rawswan_filt2.rds")
  # mvals = getM(rawswan_filt2)
  # bvals = getBeta(rawswan_filt2)


  ## Run pathways for each comparison
  for (comparison in comparisons {
    # Get names of samples used in comparison
    comparison.groups <- names(which(cont.matrix[,comparison] != 0))
    samples <- TARGETS[TARGETS$Sample_Group %in% comparison.groups, "ID"]

    # Get comparison data
    resall <- topTable(fit2, adjust="BH", num=Inf, coef=comparison)

    ## If already generated
    #resall$genename <- allgenelabels

    # Annotate all probes with geneid
    ann450k_tmp <- ann450k[row.names(resall),]
    allgenes <- strsplit(ann450k_tmp$UCSC_RefGene_Name, ";", fixed=TRUE)
    allgenelabels <- sapply(allgenes, "[", 1)
    resall$genename<-allgenelabels

    # Get entrez vector in order of external gene names
    res_entrez_ids = gene_name_to_entrez[match(resall$genename, gene_name_to_entrez$external_gene_name), "entrezgene_id"]
    anno = data.frame(res_entrez_ids)
    rownames(anno) = rownames(resall)
    colnames(anno) = "Entrez ID"

    # Get entrezID to match to reactome pathway
    resall$`Entrez ID` <- anno$`Entrez ID`
    resallf<-resall[which(resall$`Entrez ID`!="NA"),]
    res2 <- resallf[ resallf$`Entrez ID` %in% keys( reactome.db, "ENTREZID" ) & !is.na( resallf$adj.P.Val) , ]
    cat("nrows is: ")
    cat("ID's after filter ", nrow(res2), "\n")
    cat("ID's before FDR filter ", nrow(resallf),"\n")

    ## Take unique records, shoudl give same result
    #res2$entrezid <- as.character(substring(make.names(as.integer(res2$`Entrez ID`), unique=T), 2))
    #cat("Unique character res2 ", length(unique(as.character(res2$entrezid))), "\n")

    ## Use all probes(not just unique genes)
    res2$entrezid <- res2$`Entrez ID`

    ## Annotate subset of DB (make unique keys)
    reactomeTable <- AnnotationDbi::select( reactome.db, keys=as.character(res2$entrezid), keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID"))
    reactomeTable <- unique(reactomeTable)

    ## Convert Entrezt back to common gene name
    geneTable <- unique(AnnotationDbi::select(org.Hs.eg.db, keys=reactomeTable$ENTREZID, keytype="ENTREZID", columns=c("SYMBOL","ENTREZID")))
    reactomeDB_human<-merge(reactomeTable, geneTable, by="ENTREZID", all.x=T)[2:3]

    ## Make a vector of genes in each pathway
    reactomeDB_human_uniq <- unique(reactomeDB_human)
    res<-aggregate(SYMBOL~REACTOMEID, reactomeDB_human_uniq, FUN= paste, collapse='|')

    ## Final DF (matrix of reactome ID (rows) and entrezid (cols), TRUE/FALSE)
    incm <- do.call( rbind, with(reactomeTable, tapply(
      ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x ) ))

    #colnames(incm) <- res2$entrez
    colnames(incm) <- rownames(res2)

    ##remove paths with <20 genes or +80 genes
    within <- function(x, lower, upper) (x>=lower & x<=upper)
    incm <- incm[ within(rowSums(incm), lower=20, upper=80), ]

    ## Run t-test
    testCategory <- function( reactomeID ) {
      isMember <- incm[ reactomeID, ]
      possibleError<-tryCatch(t.test( res2$logFC[ isMember ] )$p.value, error=function(e) e)
      data.frame(
        reactomeID = reactomeID,
        numGenes    = sum( isMember ),
        avgLFC      = mean( res2$logFC[isMember] ),
        sdLFC       = sd( res2$logFC[isMember] ),
        zValue      = mean( res2$logFC[isMember] ) /sd( res2$logFC[isMember] ),
        strength    = sum( res2$logFC[isMember] ) / sqrt(sum(isMember)),
        pvalue=if(!inherits(possibleError, "error")){t.test( res2$logFC[ isMember ] )$p.value} else {"NA"},
        #pvalue      = t.test( res2$logFC[ isMember ] )$p.value,
        reactomeName= if(is.character(reactomePATHID2NAME[[reactomeID]])){reactomePATHID2NAME[[reactomeID]]} else {print("NA")},
        #reactomeName = reactomePATHID2NAME[[reactomeID]][1],
        stringsAsFactors = FALSE, check.rows=T )
    }

  tmp<-matrix(data = NA, ncol=8, nrow=nrow(incm)*2)
  k=1
  for (i in rownames(incm)) {
    tmp[k,] <- as.character(testCategory(i)[1,])
    k=k+1
  }

    results<-data.frame(tmp[1:nrow(incm), ])
    colnames(results) <- c("reactomeID","numGenes","avgLFC","sdLFC","zValue","strength","pvalue","reactomeName")
    colnames(res)<-c("reactomeID","symbol")
    results$padj <- p.adjust( as.numeric(as.character(results$pvalue)), "BH" )
    results<-dplyr::inner_join(results, res, by="reactomeID", all=T)
    pathway_results<-results[order(results$pvalue,decreasing = T), ]
    pathway_results_p05<-pathway_results[which(as.numeric(as.character(pathway_results$pvalue))<0.05), ]
    pathway_results_LFC<-head(results[order(abs(as.numeric(as.character(results$avgLFC))), decreasing=T), ], 50)

    write.table(pathway_results, paste(CONTRACT,"_", comparison, "_Reactome", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
    write.table(pathway_results_p05, paste(CONTRACT, "_", comparison, "_Reactome_Sig05", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
    write.table(pathway_results_LFC, paste(CONTRACT, "_", comparison, "_Reactome_top50_avgLFC", ".csv", sep=""),
              quote=F, col.names=T, row.names=F, sep = ",")
  }
}

