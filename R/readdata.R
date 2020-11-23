#' A function to load array data
#'
#' @param contractid The name of the contract
#' @param samplesheet The suffix of the sample sheet
#' @param workdir The path of the working directory
#' @return The readdata object containing the methylation experiment and accompanying annotations
#' @export
arrayread <- function(contractid, workdir, samplesheet) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package \"minfi\" needed for this function to work. Please install it.",
    call. = FALSE)
  }

  if (!requireNamespace("IlluminaHumanMethylationEPICanno.ilm10b4.hg19", quietly = TRUE)) {
    stop("EPIC database is needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(missing(contractid)){
    message("error: contractid is missing…")
    stop()
  }

  if(missing(workdir)){
    message("error: workdir is missing…")
    stop()
  }

  if(missing(samplesheet)){
    message("error: samplesheet is missing…")
    stop()
  }

  ## Create dir structure
  dir.create(paste0(workdir, "/secondary_analysis"))

  print("Loading Annotations..........")
 .annotation <- c(array = "IlluminaHumanMethylationEPICanno", annotation = "ilm10b4.hg19")
 .ann450 = minfi::getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
 .ann450k <<- .ann450[,c(1:2,12,13,18,19,22:31)]

  print("Loading Data..........")
  targets <<- tryCatch(minfi::read.metharray.sheet(base=workdir, pattern=samplesheet),
                  error=function(e) {print(paste("error: Samplesheet problem detected:", samplesheet))})
  epic <<- withCallingHandlers({minfi::read.metharray.exp(base = paste0(workdir,"/", "IDAT"))},
                               error = function(e) {print(paste("error: IDAT/ not found in ", workdir, "or could not load"))})
                  # warning = function(w) {print(paste('warning: ', w))},
                  # finally = {minfi::read.metharray.exp(base = paste0(workdir,"/", "IDAT"))})

  ## Make sample names and annotations match
  print("Formatting Data..........")
  targets$Array <- gsub("_","", x=targets$Array)

  # print(slotNames(epic))
  slot(epic, "annotation") <- c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
  sampleID <- paste(targets$Slide, targets$Array, sep="_")

  ## Check they match and rename data
  print("Confirming Samples..........")
  if (!identical(sampleID, minfi::sampleNames(epic))){
    message("error: Samplesheet does not match data..")
    message("       samples in sheet: ", nrow(targets))
    message("       samples in data: ", nrow(epic))
    stop()
  }

  ## Rename samples in array
  # targets$Sample_Group = sapply(strsplit(targets$Sample_Name, "-", ), `[`, 1)
  targets$ID = paste(targets$Sample_Group, targets$Sample_Name, sep=".")
  minfi::sampleNames(epic) = targets$ID

  # Bind to list
  value <- list(sampledata=targets, arraydata=epic)
  # attr(value, 'class') <- 'RGChannelSet'
  value
  # return(targets)
  # return(epic)
  # assign("epic", epic, envir = .GlobalEnv)
  # assign("targets", targets, envir = .GlobalEnv)
}
