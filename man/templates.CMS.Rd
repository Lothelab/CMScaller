% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{templates.CMS}
\alias{templates.CMS}
\title{consensus molecular subtype (CMS) templates}
\format{An object of class \code{data.frame} with 530 rows and 3 columns.}
\usage{
templates.CMS
}
\description{
consensus molecular subtype (CMS) templates
}
\details{
Colorectal cancer Consensus Molecular Subtypes (CMS) prediction
templates for \code{\link[CMScaller]{ntp}}. Marker genes were
identified using TCGA RNA-sequencing data. \code{templates$probe} refers to
Entrez ids.
}
\examples{
head(templates.CMS)
table(templates.CMS$class)
}
\references{
Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med [Internet]. 2015 [cited 2015 Nov 5];advance online publication. Available from: \url{http://www.nature.com/nm/journal/vaop/ncurrent/full/nm.3967.html}

Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. 2010;5:e15543.
}
\keyword{datasets}