library(testthat)
library(CMScaller)

x <- stats::rnorm(10000, 1)
set.seed(1)
expect_equal(corCosine(x, -x), matrix(-1)) # equals -1
expect_equal(corCosine(x, x), matrix(1))  # equals 1
expect_gt(mean(replicate(1000, corCosine(x,sample(x)))), 0) # expectation .5

expect_error(corCosine(1:4, 4))
expect_error(corCosine(4:1, 4))

# check that cosine and Pearson's correlation distances are equivalent
expect_equal(sqrt(1/2*(1-((cor(x,-x))))), 1)
expect_equal(as.vector(CMScaller:::simToDist(corCosine(x,-x))),1)

