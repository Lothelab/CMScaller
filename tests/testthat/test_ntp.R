library(testthat)
library(Biobase)
library(CMScaller)

# courtsey of Tyler Rinker http://stackoverflow.com/questions/24252241
is.bad <- function(code) {
    isTRUE(tryCatch(code,
                    error = function(c) TRUE,
                    warning = function(c) TRUE
    ))
}

set.seed(1)
emat <- exprs(crcTCGAsubset)

# checkInput ##################################################################

templates.phony <- templates.MSI
templates.phony$probe <- as.factor(templates.phony$symbol)
expect_true(is.bad(ntp(emat, templates.phony, nPerm=1)))

# ematAdjust
emat <- ematAdjust(emat)
expect_lte(mean(emat, na.rm = TRUE), 1e-2)
expect_lte(sd(emat, na.rm = TRUE), 1.1)

# checkOutput #################################################################

res <- ntp(emat, templates.MSI, nPerm=1)
expect_named(res)
expect_equal(rownames(res), colnames(emat))
expect_equal(rownames(res), colnames(emat))
expect_equivalent(as.numeric(res$prediction),
                apply(res[,grepl("d\\.", colnames(res))], 1, which.min))
expect_equal(levels(res$prediction), levels(templates.MSI$class))

# checkPvalues ################################################################

# no permutations
res <- ntp(emat, templates.MSI, nPerm=1)
expect_equal(min(res$p.value), 1)
# permuations
res <- ntp(emat, templates.MSI, nPerm=500)
expect_lte(quantile(res$FDR, probs = .5), .1)
# random templates
set.seed(1)
templates.random <- templates.MSI
templates.random$probe <- sample(rownames(emat), size=nrow(templates.random))
res <- ntp(emat, templates.random, nPerm=500)
expect_gte(min(res$FDR), .1)

# check seeds
res1 <- ntp(emat, templates.random, seed = 1)
res2 <- ntp(emat, templates.random, seed = 1, nCores = 4L)
all(res1$p.value==res2$p.value)

# random expression and templates
N <- 1000;P <- 5000;K <- 4;nPerm=1000
rowN <- make.names(seq_len(P))

dos <- function(n) sum(sample(1:100, size = 5))

X <- matrix(rnorm(P*N, 0, 1), ncol=N)
T <- matrix(rbinom(P*K, 1, 1e-2), ncol=K)
rownames(X) <- rowN
temp <- lapply(seq_len(K), function(k) rowN[T[,k] ==1])
names(temp) <- paste0("k", seq_len(K))
temp <- ntpMakeTemplates(temp, resDEG=FALSE)
p <- ntp(X, temp,nPerm=nPerm)$p.value
expect_warning(expect_error(CMScaller(X, templates.CMS)))
q <- seq(0,1,by=.1)
expect_true(all(round(q,1) == round(quantile(p,probs = q),1)))

# sample-template distances
# estimate prediction confidence for sample n(

# checkFunctions ##############################################################

res <- subSetNA(res, FDR = .1)
expect_true(all(is.na(res$prediction)))
res <- subSetNA(res, pValue = 1, FDR = 1)
expect_true(all(!is.na(res$prediction)))

