% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gwscaR.R
\name{infer.mat.alleles}
\alias{infer.mat.alleles}
\title{A function to infer maternal alleles from a vcf file}
\usage{
infer.mat.alleles(dad.kid, vcf)
}
\arguments{
\item{dad.kid}{A data.frame with two columns, each one with matching father and offspring IDs}

\item{vcf}{A data.frame with genomic data in vcf format}
}
\value{
A data.frame of inferred maternal alleles
}
\description{
A function to infer maternal alleles from a vcf file
}
