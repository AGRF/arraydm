#' A wrapper to run all functions in sequence
#'
#' @param contractid Character. The name of the contract or other value
#' @param contrastlist List. A named list of comparisons to undertake.
#' @param sampledata A data frame with the samplesheet.Generate manually or in arraydm::readdata().
#' @param workdir Where to save the plots. (Default: working directory)
#' @return Full suite of differential methylation analysis
#' @export
auto_arraydm <- function (contractid, sampledata, workdir, contrastlist) {
  if(missing(contractid)){
    message("error: contractid is missing…")
    stop()
  }

  if(missing(sampledata)){
    message("error: sampledata is missing…")
    stop()
  }

  if(missing(workdir)){
    message("error: workdir is missing…")
    stop()
  }

  ## Set env
  options(scipen=999)

  ## Open a file for stderrs and set error logging options
  errorfile <- file(paste0(workdir, "/errors.Rout"), open="wt")
  sink.number(type = c("output", "message"))
  sink(errorfile, append=TRUE, split=FALSE)

  # Write logfile header
  writeLines(paste0("Running DM pipeline version ", packageVersion("arraydm"), "; Date/Time: ", Sys.time()), errorfile)

  # Step logging function
  logstep <- function(step) {
    writeLines(paste0("\n", "Step: ", step, "\n------------------"), errorfile)
  }

  ## Read data - takes a few minutes
  step = "Reading inputs"
  logstep(step)
  data <- arrayread(contractid=contractid, workdir=workdir, samplesheet=".csv$")

  ## Pre-QC
  step = "Pre-QC"
  logstep(step)
  arraypreqc(contractid=contractid, arraydata=epic, sampledata=targets, workdir=workdir)

  ## Normalisation
  step="Normalisation"
  logstep(step)
  arraynorm(contractid=contractid, arraydata=epic, mdataraw=rawmeth, sampledata=targets, workdir=workdir)

  ## Post-QC
  step="Post-QC"
  logstep(step)
  arraypostqc(contractid=contractid, swandata=swannorm, pdata=pvals_filt, sampledata=targets, workdir=workdir)

  ## Variance Analysis
  step="Variance Analysis"
  logstep(step)
  arrayvar(contractid=contractid, mdata=mvals, sampledata=targets, workdir=workdir)

  ## Run Differential Methylation Tests
  step="Differential Methylation"
  logstep(step)
  arraymodel(contractid=contractid, sampledata=targets, workdir=workdir, mdata=mvals, contrastl=contrast_list)

  ## Array Summarisation
  step="Summarise Results"
  logstep(step)
  arraysumm(contractid=contractid, bdata=bvals , sampledata=targets, workdir=workdir, contrastm=cont.matrix)

  ## Pathways Analysis
  step="Pathways Analysis"
  logstep(step)
  arraypaths(contractid=contractid, sampledata=targets, workdir=workdir, myrds=lrt_list.RDS, contrastm=cont.matrix)

  ## Gene Ontology
  step="Gene Ontology"
  logstep(step)
  arraygo(contractid=contractid, sampledata=targets, workdir=workdir, myrds=lrt_list.RDS, contrastm=cont.matrix)

  # Write logfile footer
  writeLines(paste0("Finished DM pipeline at Date/Time: ", Sys.time()), errorfile)

  ## Close error logging
  # close(errorfile)
}
