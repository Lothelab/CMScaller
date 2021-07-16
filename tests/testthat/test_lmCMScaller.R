library(CMScaller)
library(testthat)

## test that classifications in paper are reproduced

data("mcrcOSLOsubset")
res <- lmCMScaller(mcrcOSLOsubset, posterior=.6)
cross_tab <- table(RF=res$prediction, PAM=mcrcOSLOsubset$cms_pam, useNA="alw")
cross_tab
test_that('tests that predictions are better than expectation', {
    expect_gte(CMScaller:::accuracy(res$prediction, mcrcOSLOsubset$cms_pam), .95)
})

