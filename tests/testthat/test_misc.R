library(testthat)
library(CMScaller)

classCol <- getOption("subClassCol")
expect_equal(class(classCol), "character")

# internal functions
expect_equal(CMScaller:::distToSim(CMScaller:::simToDist(.5)), .5)
expect_true(CMScaller:::packageExists("CMScaller"))

# documentation
expect_silent(citation("CMScaller"))
expect_silent(news(package = "CMScaller"))

