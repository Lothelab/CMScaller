% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ntpMakeTemplates.R
\name{ntpMakeTemplates}
\alias{ntpMakeTemplates}
\title{nearest template prediction templates}
\usage{
ntpMakeTemplates(
  geneSets,
  resDEG = TRUE,
  lfc = 1,
  padj = 0.1,
  topN = NULL,
  verbose = getOption("verbose", FALSE)
)
}
\arguments{
\item{geneSets}{list of vectors with gene identifiers (coerced to
character).}

\item{resDEG}{logical, if \code{geneSets} is result from \link{subDEG}.}

\item{lfc}{numeric, log2fold-change threshold (only applicable if
\code{resDEG=TRUE}).}

\item{padj}{numeric, \eqn{p_{adj}}-value threshold (only applicable if
\code{resDEG=TRUE}).}

\item{topN}{integer, (maximum) number of genes per class.}

\item{verbose}{logical, whether console messages are to be displayed.}
}
\value{
A data.frame formated as templates for \code{\link{ntp}}.
}
\description{
Reformats list of gene sets or results from
\code{\link{subDEG}} into prediction templates for \code{\link{ntp}}.
}
\examples{
deg <- subDEG(crcTCGAsubset, crcTCGAsubset$CMS, doVoom=TRUE, sortBy="B")
tmp <- ntpMakeTemplates(deg, topN=25)
table(tmp$class)
}
\seealso{
\code{\link{ntp}}, \code{\link{subDEG}}
}
