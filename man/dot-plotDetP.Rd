% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotdetp.R
\name{.plotDetP}
\alias{.plotDetP}
\title{Visualises the detection P value of array prior to filtering}
\usage{
.plotDetP(contractid, pvals, sampledata, workdir)
}
\arguments{
\item{contractid}{Character. The name of the contract}

\item{pvals}{Data Frame. Detection P-values}

\item{sampledata}{Data Frame. The project samplesheet.Generate manually or in arraydm::readdata()}

\item{workdir}{Character. Path to output location. (Default is working directory)}
}
\value{
A barplot of the detection P values for the unfiltered array
}
\description{
Visualises the detection P value of array prior to filtering
}
