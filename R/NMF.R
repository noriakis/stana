
#' NMF
#'
#' decompose SNV, gene content, or gene family abundance to
#' strain x sample and strain x feature matrix.
#'
#' @param stana stana object
#' @param species candidate species ID
#' @param rank rank of NMF
#' @param target KO, gene or snv, default to KO
#' @param seed random seed
#' @param method NMF method, default to snmf/r
#' @param plotHeatmap plot the heatmap of presence of strain across sample
#' @import NMF
#' @export
NMF <- function(stana, species, rank=5, target="KO", seed=53, method="snmf/r",
    deleteZeroDepth=TRUE, beta=0.01, plotHeatmap=TRUE) {
    if (length(species)>1) {stop("NMF accepts only one species per run")}
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
    cat("Original features:", dim(mat)[1], "\n")
    cat("Original samples:", dim(mat)[2], "\n")
    mat <- mat[apply(mat, 1, function(x) sum(x)!=0),]
    mat <- mat[,apply(mat, 2, function(x) sum(x)!=0)]

    cat("Filtered features:", dim(mat)[1], "\n")
    cat("Filtered samples:", dim(mat)[2], "\n")

    ## Test multiple ranks
    cat("Rank", rank, "\n")
    if (method %in% c("snmf/l", "snmf/r")) {
        res <- NMF::nmf(mat, rank = rank, seed = seed, method=method, beta=beta)
    } else {
        res <- NMF::nmf(mat, rank = rank, seed = seed, method=method)
    }
    coefMat <- coef(res)
    ## Plot by default
    relab <- apply(coef(res), 2, function(x) x / sum(x))
    cat("Mean relative abundances:", apply(relab, 1, mean), "\n")

    basisMat <- basis(res)
    cat("Present feature per strain:", apply(basisMat!=0, 2, function(x) sum(x)), "\n")
    stana@NMF[[species]] <- res
    return(stana)
}

#' alphaDiversityWithinSpecies
#' @export
alphaDiversityWithinSpecies <- function(stana, species, method="shannon") {
    if (is.null(stana@NMF[[species]])) {
        stana <- NMF(stana, species)
    }
    res <- stana@NMF[[species]]
    W <- coef(res)
    if (method=="spc") {
        div <- apply(coef(stana@NMF[[1]])==0, 2, function(x) sum(x))
    } else {
        div <- vegan::diversity(t(W), index=method)
    }
    
    if (!is.null(stana@cl)) {
        nm <- listToNV(cl)
        div <- data.frame(div)
        colnames(div) <- "alpha_diversity"
        div[["group"]] <- nm[row.names(div)]
    }
    return(div)
}