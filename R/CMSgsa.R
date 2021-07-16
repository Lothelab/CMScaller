#' CMS Gene Set Analysis
#' @export
#' @description \code{\link[limma]{camera}} gene set analysis (GSA) against 15
#' pre-selected CMS-informative gene sets
#' (\code{\link{geneSets.CMS}}).
#' @export
#' @param emat a numeric expression matrix with sample columns, gene rows and
#' Entrez rownames. Microarray data should be normalized and log2-transformed.
#' For RNA-seq data, raw counts or RSEM values could be used directly by setting
#' \code{RNAseq=TRUE}.
#' @param class a factor vector specifying sample classes.
#' @param RNAseq a logical, set to TRUE if emat is untransformed, non-normalized
#' sequencing counts or RSEM values.
#' @param ... additional arguments passed to \code{\link{subCamera}}.
#' @details See \code{\link{subCamera}} for output details.
#' @return a heatmap and \code{\link[limma]{camera}} output (list,
#' invisible). In heatmap, red and blue indicates relative up- and
#' down-regulation respectively. Color saturation reflects significance.
#' Nominal `camera` \eqn{p}-values are used as input for visualization.
#' @seealso \code{\link[limma]{camera}}, \code{\link{subCamera}}, \code{\link{geneSets.CMS}}
#' @references Guinney J, Dienstmann R, Wang X, de Reyniès A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350–6. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4636487/}
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucl. Acids Res. 2015;gkv007. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/}
#' @references Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26:139-40. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/19910308/}
#' @references Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology. 2014;15:R29. Available from: \url{https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/24485249/}
#' @examples
#' cam <- CMSgsa(emat=crcTCGAsubset, class=crcTCGAsubset$CMS, RNAseq=TRUE)
#' lapply(cam, head)
CMSgsa <- function(emat, class, RNAseq=FALSE, ...) {
  if(length(class)!=ncol(emat))
    stop ("length(class)) not equal to ncol(emat))")
  subCamera(emat, class, topN=15, doVoom=RNAseq, ...,
            rowCluster = FALSE, geneList = CMScaller::geneSets.CMS)
}
