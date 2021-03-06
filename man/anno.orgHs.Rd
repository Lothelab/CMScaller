% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{anno.orgHs}
\alias{anno.orgHs}
\title{annotation for translating between gene identifiers}
\format{
An object of class \code{data.frame} with 60609 rows and 7 columns.
}
\usage{
anno.orgHs
}
\description{
Table for translating between
\href{Ensembl GENCODE}{https://www.gencodegenes.org/human/},
\href{NCBI Entrez}{https://www.ncbi.nlm.nih.gov/gene}, and
\href{HGNC symbol}{https://www.genenames.org/} human gene identifiers. Gene
biotype, exonic gene length and proportion GC-bases are also included. For
more exhaustive annotation including gene names, GO terms and so on, please
refer to R Bioconductor packages \code{AnnotationDbi} and \code{org.Hs.eg.db}.
}
\details{
Database versions are indicated in package startup message. NCBI
Entrez identifiers matched using HGNC symbol against
\code{\link[org.Hs.eg.db]{org.Hs.eg.db}} package.
}
\examples{
head(anno.orgHs)
dim(anno.orgHs)
colSums(apply(anno.orgHs, 2, is.na))
colSums(!apply(anno.orgHs[,1:4], 2, duplicated))
# some extremes
anno.orgHs[which.max(anno.orgHs$length),]
anno.orgHs[which.min(anno.orgHs$length),]
anno.orgHs[which.max(anno.orgHs$GC),]
anno.orgHs[which.min(anno.orgHs$GC),]
attributes(anno.orgHs)$genome
}
\references{
Braschi B, Denny P, Gray K, Jones T, Seal R, Tweedie S, et al. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. 2019;47:D786–92. Available from \url{https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/30304474/}

Carlson M. org.Hs.eg.db: Genome wide annotation for Human. R package version 3.11.4. 2020. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3159387/}

Frankish A, Diekhans M, Ferreira A-M, Johnson R, Jungreis I, Loveland J, et al. GENCODE reference annotation for the human and mouse genomes. Nucleic Acids Res. 2019;47:D766–73. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6323946/}

Maglott D, Ostell J, Pruitt KD, Tatusova T. Entrez Gene: gene-centered information at NCBI. Nucleic Acids Res. 2007;35:D26-31. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/17148475/}
}
\seealso{
\code{\link{fromTo}}, \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
}
\keyword{datasets}
