#' liver metastases CMS classification
#' @export
#' @param emat a numeric expression matrix with sample columns, gene rows and
#' Entrez rownames. Microarray data should be normalized. For RNA-seq data,
#' counts or RSEM values could be used directly by setting \code{RNAseq=TRUE}.
#' @param rowNames a character, either "entrez" (NCBI Entrez),
#' "symbol" (HGNC symbol) or "ensg" (Ensembl). If set to other than "ensg",
#' \code{\link{replaceGeneId}} is used to translate \code{rownames(emat)}.
#' @param RNAseq a logical, set to TRUE if emat is untransformed, non-normalized
#' sequencing counts or RSEM values.
#' @param posterior numeric, sets minimum prediction confidence threshold.
#' @details \code{lmCMScaller} provides CMS classification of colorectal cancer
#' liver metastasis based on a pre-defined \code{\link[randomForest]{randomForest}}
#' model. If \code{RNA-seq=TRUE}, a pseudocount of 0.25 is added, matrix
#' log2-transformed and quantile normalized (\code{\link[limma]{normalizeQuantiles}})
#' prior to scaling/centering and prediction.
#' @return a data frame with columns of class-wise posterior probabilities,
#' predictions and maximum posterior probilities. Rownames equal \code{emat}
#' colnames.
#' @references Guinney J, Dienstmann R, Wang X, de Reyniès A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350–6. Available from: \url{https://www.nature.com/articles/nm.3967}
#' @references Liaw A, Wiener M. Classification and Regression by randomForest. R News. 2002;2:18–22. Available from: \url{https://cran.r-project.org/doc/Rnews/Rnews_2002-3.pdf}
#' @examples
#' library(randomForest)
#' data("mcrcOSLOsubset")
#' emat <- mcrcOSLOsubset[,!duplicated(mcrcOSLOsubset$`Patient ID`)]
#' res1 <- lmCMScaller(emat)
#' sum(is.na(res1))
#' head(res1)
#' table(res1$prediction, emat$cms_pam)
#' ## see that it handles (some) missing values
#' res2 <- lmCMScaller(emat[-c(1:100),])
#' sum(is.na(res2))
#' table(res1$prediction, res2$prediction)
lmCMScaller <- function(emat, rowNames="entrez", RNAseq=FALSE, posterior=0.5) {

    # checkInput ##############################################################

    # check datatype input and try to coerce to matrix
    if (class(emat)[1] == "ExpressionSet") {
        emat <- suppressPackageStartupMessages(Biobase::exprs(emat))
    }
    if (class(emat)[1] == "data.frame") emat <- as.matrix(emat)
    if (is.vector(emat)) emat <- matrix(emat, dimnames = list())
    if (is.null(rownames(emat))) stop("missing Ensembl id rownames(emat)")

    if (ncol(emat) < 100) warnings("too few samples - high prediction variance",
                                  call.=FALSE)

    if (rowNames != "entrez") {
        if (!rowNames %in% c("symbol", "ensg"))
            stop("invalid rowNames, must be either entrez, symbol or ensg")
        emat <- replaceGeneId(emat, id.in=rowNames, id.out="entrez")
    }

    # log2-transform and quantile normalize RNA-seq data
    if (isTRUE(RNAseq)) {
            message("performing log2-transform and quantile normalization...")
        emat <- limma::normalizeQuantiles(log2(emat+.25))
    }

    # sanity check - whether rownames appear to be Entrez ids
    is.na.rows <- is.na(fromTo(rownames(emat), rough=TRUE))
    mm <- sum(is.na.rows)/nrow(emat)
    if (mm > 0.15) {
        message (paste0(sum(is.na.rows),"/",nrow(emat),
                        " rownames(emat) failed to match to human gene identifiers"))
        warning (paste0("verify that rownames(emat) are ", rowNames),
                 call.=FALSE)
    }

    ## NA impute missing features
    features <- CMScaller:::subData$model_mCMS$xNames
    missing <- setdiff(features, rownames(emat))

    if(length(missing)>0) {
        message(paste(missing, collapse = "; "))
        message(paste(length(missing), 'missing features, NAs imputed'))
        ndata <- matrix(data=NA, ncol=ncol(emat), nrow=length(missing),
                        dimnames = list(missing, colnames(emat)))
        emat <- rbind(emat, ndata)
    }
    ##TODO more robust imputation - rough fix that seems to work with gene exp data
    emat <- randomForest::na.roughfix(emat)

    # scale and center data, basically a wrapper for scale() function
    emat <- t(ematAdjust(emat))

    # rfPredict ###############################################################

    pred <- randomForest:::predict.randomForest(subData$model_mCMS, emat, "prob")
    res <- data.frame(pred,
        prediction=apply(pred, 1, which.max),
        posterior=apply(pred, 1, max))
    setNA <- res$posterior < posterior
    res$prediction <- paste0("CMS", 1:4)[res$prediction]
    res$prediction <- factor(res$prediction, levels=paste0("CMS", 1:4))
    res$prediction[setNA] <- NA

    # output ##################################################################
    return(res)
}
