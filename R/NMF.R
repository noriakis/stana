
#' NMF
#' @param stana stana object
#' @param species candidate species ID
#' @param rank rank of NMF
#' @param target KO, gene or snv, default to KO
#' @param seed random seed
#' @param method NMF method, default to snmf/r
#' @import NMF
#' @export
NMF <- function(stana, species, rank=10, target="KO", seed=53, method="snmf/r",
    deleteZeroDepth=TRUE, beta=0.01) {
    if (target=="KO") {
        mat <- stana@kos[[species]]
    } else if (target=="genes") {
        mat <- stana@genes[[species]]
    } else if (target=="snps") {
        mat <- stana@snps[[species]]
        if (deleteZeroDepth) {
           mat <- mat[rowSums(mat == -1)==0,]
           qqcat("After filtering `-1`, position numbers: @{dim(mat)[1]}\n")
        }
    }
    res <- NMF::nmf(mat, rank = rank, seed = seed, method=method, beta=beta)
    stana@NMF[[species]] <- res
    return(stana)
}

#' alphaDiversityWithinSpecies
#' @export
alphaDiversityWithinSpecies <- function(stana, species) {
    if (is.null(stana@NMF[[species]])) {
        stana <- NMF(stana, species)
    }
    res <- stana@NMF[[species]]
    W <- coef(res)
    div <- vegan::diversity(t(W))
    
    if (!is.null(stana@cl)) {
        nm <- listToNV(cl)
        div <- data.frame(div)
        colnames(div) <- "alpha_diversity"
        div[["group"]] <- nm[row.names(div)]
    }
    return(div)
}