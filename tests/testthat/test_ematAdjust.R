library(testthat)
library(Biobase)
library(CMScaller)


emat <- exprs(crcTCGAsubset)[1:100,]
emat <- ematAdjust(emat, normMethod = "quantile")

expect_equal(sd(emat[1,]), 1)
expect_equal(mean(emat[1,]), 0)
expect_lte(max(emat), 10)

emat <- exprs(crcTCGAsubset)[1:100,]
emat <- ematAdjust(emat, normMethod = "quantile", scale=FALSE, center=TRUE)
expect_equal(mean(emat[1,]), 0)
emat <- ematAdjust(emat, scale=TRUE, center=TRUE)
expect_equal(sd(emat[1,]), 1)
