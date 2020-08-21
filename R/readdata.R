#' A function to load array data
#'
#' @param contractid The name of the contract
#' @param workdir The path of the working directory
#' @param samplesheet The suffix of the sample sheet
#' @return The readdata object containing the methylation experiment
#' @export
readdata <- function(contractid, workdir, samplesheet) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
    call. = FALSE)
  }

  if(missing(contractid)){
    message("contractid is missing…")
  }

  if(missing(workdir)){
    message("workdir is missing…")
  }

  if(missing(samplesheet)){
    message("samplesheet is missing…")
  }

  targets <- tryCatch(read.metharray.sheet(base=workdir, pattern=samplesheet),
                  error=function(e) {print(paste("Samplesheet error in", samplesheet))})
  data <- tryCatch(read.metharray.exp(base = paste0(workdir,"/", IDAT"),
                  error=function(e) {print(paste("IDAT/ not found in ", workdir))})
  value <- list(targets=targets, data=data)
  attr(value, 'class') <- 'rawdata'
  value
}

## Usage:
# readdata(contractid="ILMLEPIC-15325", workdir="/data/Analysis/Genotyping/CAGRF12345/", samplesheet=".csv$")




