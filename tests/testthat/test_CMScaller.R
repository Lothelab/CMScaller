library(testthat)
library(Biobase)
library(CMScaller)

# checkPredictions ############################################################

# provides checking of internal voomTransform(), templates.CMS and ntp()
crcTCGAsubset <- crcTCGAsubset[,!crcTCGAsubset$isTrain & crcTCGAsubset$platform == "hiseq"]
templates.CMS <- templates.CMS[!duplicated(templates.CMS$probe),]
set.seed(1)
res <- CMScaller(crcTCGAsubset, RNAseq=TRUE, doPlot = FALSE)
expect_gte(CMScaller:::accuracy(res$prediction, crcTCGAsubset$CMS), .83)

# check geneTranslations ######################################################
emat.symbol <- replaceGeneId(crcTCGAsubset, id.in="entrez", id.out = "symbol")
emat.entr <- replaceGeneId(crcTCGAsubset, id.in="entrez", id.out = "entrez")
emat.ensg <- replaceGeneId(crcTCGAsubset, id.in="entrez", id.out = "ensg")
expect_true(all(emat.symbol==emat.entr))

expect_error(CMScaller(emat.ensg, rowNames = "ncbi", RNAseq = TRUE))
res.symb <- CMScaller(emat.symbol, rowNames = "symbol", RNAseq = TRUE)
res.entr <- CMScaller(emat.entr, rowNames = "entrez", RNAseq = TRUE)
expect_equal(CMScaller:::accuracy(res.symb$pred, res$prediction), 1)
expect_equal(CMScaller:::accuracy(res.symb$pred, res.entr$prediction),1)

