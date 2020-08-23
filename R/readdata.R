#' A function to load array data
#'
#' @import arraydm
#' @param contractid The name of the contract
#' @param workdir The path of the working directory
#' @param samplesheet The suffix of the sample sheet
#' @return The readdata object containing the methylation experiment
#' @export
readdata <- function(contractid, workdir, samplesheet) {
  if (!requireNamespace("minfi", quietly = TRUE)) {
    stop("Package \"minfi\" needed for this function to work. Please install it.",
    call. = FALSE)
  }

  if(missing(contractid)){
    message("contractid is missing…")
    exit()
  }

  if(missing(workdir)){
    message("workdir is missing…")
    exit()
  }

  if(missing(samplesheet)){
    message("samplesheet is missing…")
    exit()
  }

  targets <- tryCatch(minfi::read.metharray.sheet(base=workdir, pattern=samplesheet),
                  error=function(e) {print(paste("error: Samplesheet problem detected:", samplesheet))})
  data <- tryCatch(minfi::read.metharray.exp(base = paste0(workdir,"/", "IDAT")),
                   error = function(e) {print(paste("error: IDAT/ not found in ", workdir, "or could not load"))},
                   warning = function(w) {print(paste('warning: ', w));  })
  value <- list(targets=targets, data=data)
  attr(value, 'class') <- 'rawdata'
  value
}
