genexp <- function(eset, symbol="MIR31HG") {
    return(c(Biobase::exprs(eset[
        fromTo(symbol,id.in="symbol", id.out="ensg"),])))
}
