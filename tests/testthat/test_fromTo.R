library(testthat)
library(CMScaller)

# fromTo
p <- seq_len(10)
m <- anno.orgHs
x <- anno.orgHs$entrez[p]

# expect true
expect_true(all(fromTo(x)==m$symbol[p]))
expect_true(all(x == fromTo(fromTo(x))))
expect_is(x, "character")

# expect error
expect_error(fromTo(NA))
expect_error(fromTo(c(NA, 2, 3)))
expect_error(fromTo(9E9))
