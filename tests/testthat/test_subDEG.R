library(testthat)
library(Biobase)
library(CMScaller)

set.seed(1)
emat <- crcTCGAsubset
class <- crcTCGAsubset$CMS.Syn
deg <- subDEG(emat, class, doVoom = TRUE)

# checkOutput #################################################################
expect_equal(names(deg), levels(class))
cms4 <- deg$CMS4
markers <- c("FERMT2", "VCAN", "LUM")
expect_true(all(markers %in% head(cms4$symbol[cms4$logFC>0], 250)))
