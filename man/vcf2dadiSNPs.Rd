% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf2dadi.R
\name{vcf2dadiSNPs}
\alias{vcf2dadiSNPs}
\title{Convert a vcf file to a SNPs input file for dadi}
\usage{
vcf2dadiSNPs(vcf, filename = NULL, pop.list = NA, pop.map = NULL,
  whitelist = NULL, outgroup = NA)
}
\arguments{
\item{vcf}{The vcf dataframe (from parse.vcf) or the name of the vcf file}

\item{filename}{The name fo the file you want the data written to. If not provided, no file is written}

\item{pop.list}{Optional list of population names. If this is provided, individuals should have the population names in their IDs. If neither pop.map or pop.list are provided, all individuals are placed in one population.}

\item{pop.map}{Optional map of individual IDs in column 1 and their population names in column 2 (same as a Stacks population map). Can either be a file or an R object. If neither pop.map or pop.list are provided, all individuals are placed in one population.}

\item{whitelist}{Optional list of SNPs to include in the output file. If not provided, all SNPs in the vcf will be output.}

\item{outgroup}{Optional 2-column dataframe containing ingroup (column 1) and outgroup (column 2) sequences. If not provided, both will be "---".}
}
\value{
a dataframe containing the contents of the vcf file, including headers.
}
\description{
Convert a vcf file to a SNPs input file for dadi
}
\examples{
vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
dadi<-vcf2dadiSNPs(vcf.file)
}
\seealso{
Gutenkunst et al. 2013
}
