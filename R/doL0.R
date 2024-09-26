#' doL0
#' 
#' Feature selection of classification between category
#' using snv frequency and gene copy number matrix, based on L0Learn
#' 
#' @param stana stana object
#' @param species candidate species ID
#' @param cl named list of category
#' @param target default to snps
#' @param mat if target is not snps, provide preprocessed gene matrix
#' otherwise the raw gene matrix is used.
#' @param deleteZeroDepth delete zero depth SNV
#' @param argList passed to L0Learn()
#' @return L0Learn object
#' @import ggplot2
#' @export
doL0 <- function(stana, species, cl=NULL,
    target="genes", mat=NULL, argList=list(),
    deleteZeroDepth=FALSE) {
    if (length(argList)==0) {
    	qqcat("Penalty is set to L0L2 by default")
    	argList[["penalty"]] <- "L0L2"}
    if (is.null(cl)) {
        qqcat("Using grouping from the slot\n")
        cl <- stana@cl
    }
    if (target=="genes" & is.null(mat)){
        qqcat("If needed, please provide preprocessed matrix of genes to `mat`\n")
        filtDf <- stana@genes[[species]]
    } else if (target=="ko") {
    	filtDf <- stana@kos[[species]]
    } else if (target=="snps"){
        filtDf <- stana@snps[[species]]
        if (deleteZeroDepth) {
            filtDf <- filtDf[rowSums(filtDf==-1)==0,]
            qqcat("After filtering `-1`, position numbers: @{dim(filtDf)[1]}\n")
        }
    } else {
        qqcat("Proceeding with provided matrix\n")
        filtDf <- mat
    }
    qqcat("Feature number: @{dim(filtDf)[1]}\n")
    transDf <- data.frame(t(filtDf), check.names=FALSE)
    transDf <- transDf[intersect(row.names(transDf), cl |> unlist() |> unique()),]
    gr <- NULL
    for (cn in rownames(transDf)){
        for (clm in seq_along(cl)){
            if (cn %in% cl[[clm]]) {
                gr <- c(gr, names(cl)[clm])
            }
        }
    }
  
    
    qqcat("Performing L0Learn\n")
    argList[["x"]] <- transDf |> as.matrix()
    argList[["y"]] <- as.factor(gr)
    l0res <- do.call("L0Learn.cvfit", argList)
    return(l0res)
}
