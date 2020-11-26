#' Array QC
#'
#' Functions to qc array data and summarise qc results.
#'
#' @param arraydata An RGChannelSet generated manually, or in arraydm::readdata()
#' @param contractid Character. The name of the contract
#' @param sampledata Data Frame. The project samplesheet, generate manually or in arraydm::readdata()
#' @param workdir Character. Path to output location. (Default is working directory)
#' @return Plots for detection p-value and a quality filtered probe set.
#' @export
arraypreqc <- function(contractid, arraydata, sampledata, workdir=NULL) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
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

  if(missing(arraydata)){
    message("error: arraydata is missing......")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing......")
    stop()
  }

  ## Set output location
  dir.create(paste0(workdir,"/secondary_analysis/Results/qc"), showWarnings = FALSE, recursive=TRUE)

  print("Running PreQC......")

  ## Add samplenames
  sampledata$ID = paste(sampledata$Sample_Group, sampledata$Sample_Name, sep=".")
  minfi::sampleNames(arraydata) = sampledata$ID

  ## Calculate pvals
  pvals = minfi::detectionP(arraydata)

  ngroups=ifelse(length(unique(sampledata$Sample_Group)) < 3, 3, length(unique(sampledata$Sample_Group)))

  pal <<- RColorBrewer::brewer.pal(ngroups,"Set2")    # Bad practice but im doing it anyway
  print("Plotting Detection P..........")
  .plotDetP(contractid, pvals, workdir, sampledata=targets)

  ## Filter by mean detection across all samples
  samplesin = colMeans(pvals) < 0.05

  epic <- arraydata[ , samplesin]
  targets <- sampledata[samplesin, ]
  pvals <- pvals[ ,samplesin]

  n = sum(samplesin==FALSE)
  print(paste0("Dropping ", n, " samples due to detection Pvals"))

  ## Convert raw to M values and B values
  print("Generating M and B values and removing failed probes ..........")
  rawmeth <<- preprocessRaw(epic)
  rawmeth.mvals <- getM(rawmeth)
  nprobes=dim(rawmeth)[1]
  print(paste0("Number of raw probes is ", nprobes))

  print("Saving 'rawmeth' Object..........")
  saveRDS(rawmeth, paste0(workdir, "/secondary_analysis/", contractid, "_rawmeth.RDS"))

  ## Filter failed probes (null/infinite calls)
  is.na(rawmeth.mvals) <- do.call(cbind,lapply(rawmeth.mvals, is.infinite))
  rawmeth_mvals <- na.omit(rawmeth.mvals)
  nprobes=dim(na.omit(rawmeth_mvals))[1]
  print(paste0("Number of passed probes is ", nprobes))

  ## Filter Pvals
  pvals_filt <- pvals[rownames(rawmeth_mvals),]
  colnames(pvals_filt) <- colnames(rawmeth_mvals)
  identical(rownames(pvals_filt), rownames(rawmeth_mvals))

  assign("epic", epic, envir = .GlobalEnv)
  assign("targets", sampledata, envir = .GlobalEnv)
  assign("pvals", pvals, envir = .GlobalEnv)
  assign("pvals_filt", pvals_filt, envir = .GlobalEnv)
  assign("rawmeth_mvals", rawmeth_mvals, envir = .GlobalEnv)

  print("Plotting PreQC Figures..........")
  .plotLumiQC(contractid=contractid, mdata=rawmeth_mvals, sampledata=targets, pdata=pvals_filt, qctype="preQC", workdir=workdir)
}

#' A function to qc array data2
#'
#' @inheritParams arraypreqc
#' @param swandata Data Frame. The project samplesheet.Generate manually or in arraydm::readdata()
#' @param pdata Data Frame. Detection P-values.
#' @return M and B values, the RDS file and postQC plots.
#' @export
arraypostqc <- function(contractid, swandata, pdata, sampledata, workdir=NULL) {
  print("Running PostQC......")

  if (!requireNamespace("minfi", quietly = TRUE)) {
  stop("The minfi package is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    stop("EPIC database is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(contractid)){
    message("error: contractid is missing......")
    stop()
  }

  outprefix=paste0(workdir,"/secondary_analysis/Results/", contractid)

  ## set defaults
  minp = 0.01

  ## Filter by detP
  swanpvals = pdata[match(minfi::featureNames(swandata),rownames(pdata)),]
  keep = rowSums(swanpvals < minp) == ncol(swandata)
  rawswan_pfilt = swandata[!is.na(keep),]

  ## Remove known SNPs
  rawswan_filt1 = minfi::dropLociWithSnps(rawswan_pfilt)
  nprobes=(dim(rawswan_pfilt)[1])-(dim(rawswan_filt1)[1])
  print(paste0("Dropping ", nprobes, " SNP......"))

  ## Drop X, Y
  keep = !(featureNames(rawswan_filt1) %in% rownames(.ann450k)[.ann450k$chr %in% c("chrX","chrY")])
  nprobes=as.numeric(table(keep)[1])
  print(paste0("Dropping ", nprobes, " X,Y probes ......"))
  rawswan_filt2 = rawswan_filt1[keep,]

  ## Filter Pvals
  pvals_filt <- pdata[rownames(rawswan_filt2),]
  colnames(pvals_filt) <- colnames(rawswan_filt2)
  identical(rownames(pvals_filt), rownames(rawswan_filt2))
  pvals <<- pvals_filt

  ## Get M and B values
  mvals <<- minfi::getM(rawswan_filt2)
  bvals <<- minfi::getBeta(rawswan_filt2)
  #nprobes=nrow(mvals)
  #print("Keeping ", nprobes, " probes after QC")

  colnames(swanpvals) <- colnames(mvals)
  pvals_filtnorm <- swanpvals[rownames(mvals), ]

  ## Saving
  print(paste0("Saving data......"))
  saveRDS(rawswan_filt2, paste0(workdir,"/secondary_analysis/", contractid, "_rawswan_filt2.RDS"))
  write.table(mvals, paste(outprefix,"_Mvalues.csv", sep=""), quote=F, col.names=T, row.names=T)
  write.table(bvals, paste(outprefix,"_Bvalues.csv", sep=""), quote=F, col.names=T, row.names=T)

  ## Plotting post-qc
 .plotLumiQC(contractid=contractid, mdata=mvals, sampledata=targets, pdata=pvals_filtnorm, qctype="postQC", workdir)
}
