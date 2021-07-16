#' #' @export
#' #' @seealso \code{\link[pamr]{pamr.train}}
#' #' @references Tibshirani et al. (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. PNAS 99(10):6567-6572 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC124443/
#' pam <- function(X, class, threshold) {
#'     ###########################################################################
#'     ## PW Eide (2019) peteid@rr-research.no
#'     ## Lothe-group / Institute for cancer research / Oslo University Hospital
#'     ##
#'     ## implementation of pamr as described in
#'     ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC124443/
#'     ###########################################################################
#'
#'     ###########################################################################
#'     ## constants
#'     p <- nrow(X)
#'     n <- ncol(X)
#'     k <- nlevels(factor(class))
#'     n_k <- table(class)
#'     pi_k <- prop.table(n_k)
#'     C <- model.matrix(~0+class)
#'
#'     ###########################################################################
#'     ## estimate parameters
#'     xbar <- drop(X %*% rep(1,n) %*% 1/n) # overall means
#'     xbar_k <- X %*% C %*% diag(1/n_k) # class-wise means
#'     delta.xbar_k <- X - xbar_k %*% t(C) # numerator for Eq.1
#'     s <- sqrt( (delta.xbar_k^2) %*% rep(1/(n - k), n)) # Eq.2
#'     ss0 <- s+median(s)
#'     ## in the paper this is given as sqrt(1/nK-1/N)
#'     ## in pamr [v1.55] this is coded as sqrt(1/nK+1/N)
#'     m_k <- sqrt(1/n_k-1/n) # scales wihin class standard errors
#'     d_k <- (xbar_k-xbar) / (ss0 %*% (m_k)) # Eq.1
#'
#'     ###########################################################################
#'     ## make shrunken centroids
#'     d <- abs(d_k)-threshold
#'     shrunk.d_k <- sign(d_k) * d * as.numeric(d>0)
#'     shrunk.d_k <- shrunk.d_k %*% diag(m_k)
#'     non0 <- (abs(shrunk.d_k) %*% rep(1, k)) > 0
#'     if (sum(non0)<=10) {
#'         stop(paste0(sum(non0), " non-zero features"))
#'     } else {
#'         message(paste0(sum(non0), " non-zero features"))
#'     }
#'
#'     ###########################################################################
#'     ## classification function
#'     disc.pam <- function(x) {
#'         dd1 <- t( (x[non0]-xbar[non0])/ss0[non0] ) %*% shrunk.d_k[non0,]
#'         dd0 <- drop((rep(1, sum(non0)) %*%
#'                         shrunk.d_k[non0,]^2))/2-log(pi_k)
#'         drop(dd1) - dd0
#'     }
#'
#'     dd <- apply(X, 2, disc.pam)
#'     sexp <- function(x) exp(sign(x) * pmin(abs(x), .Machine$sizeof.pointer))
#'
#'     res <- data.frame(t(dd))
#'     colnames(res) <- paste0("d.", rownames(dd))
#'     res$prediction <- factor(levels(class)[apply(res, 1, which.max)],
#'                              levels=levels(class))
#'     ##
#'     res$p.value <- sapply(seq_len(n), function(i)
#'                         min(sexp(dd[,i])/drop(sexp(dd[,i]) %*% rep(1,k))))
#'     res$FDR <- p.adjust(res[,"p.value"], method="fdr")
#'
#'     return(res)
#' }
