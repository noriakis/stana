
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
alphaDiversityWithinSpecies <- function(stana, species, method="shannon", rank=5) {
    if (is.null(stana@NMF[[species]])) {
        stana <- NMF(stana, species, rank=rank)
    }
    res <- stana@NMF[[species]]
    W <- coef(res)
    if (method=="spc") {
        div <- apply(coef(stana@NMF[[1]])==0, 2, function(x) sum(x))
    } else {
        div <- vegan::diversity(t(W), index=method)
    }
    
    if (!is.null(stana@cl)) {
        nm <- listToNV(stana@cl)
        div <- data.frame(div)
        colnames(div) <- "alpha_diversity"
        div[["group"]] <- nm[row.names(div)]
    }
    return(div)
}

#' plotAbundanceWithinSpecies
#' 
#' plot abundances using factor to sample matrix produced by NMF.
#' 
#' @param stana stana object
#' @param species species ID
#' @param tss perform total sum scaling
#' @param return_data return only the data, not plot
#' @export
plotAbundanceWithinSpecies <- function(stana, species, tss=TRUE, return_data=FALSE) {
    if (is.null(stana@NMF[[species]])) {
        stana <- NMF(stana, species)
    }
    res <- stana@NMF[[species]]
    W <- coef(res)
    if (tss) {
        W <- apply(W, 2, function(x) x / sum(x))
    }
    W <- data.frame(t(W))
    if (!is.null(stana@cl)) {
        nm <- listToNV(stana@cl)
        W[["group"]] <- nm[row.names(W)]
    }
    colnames(W) <- c(as.character(seq_len(ncol(W)-1)),"group")

    if (return_data) {
        return(W)
    }
    W %>% tidyr::pivot_longer(1:(ncol(W)-1)) %>%
        ggplot(aes(x=group, y=value))+
        geom_boxplot()+
        facet_wrap(.~name)
}


#' pathwayWithFactor
#' convert KO matrix per factor obtained by NMF function
#' to pathway to factor matrix by summing the KOs in the pathway.
#' @param stana stana boject
#' @param species species ID
#' @param tss perform total sum scaling to the resulting matrix
#' @export
pathwayWithFactor <- function(stana, species, tss=FALSE) {
  dat <- stana@NMF[[species]]
  dat <- basis(dat)

  bfc <- BiocFileCache()
  url <- bfcrpath(bfc,"https://rest.kegg.jp/link/ko/pathway")
  url2 <- bfcrpath(bfc,"https://rest.kegg.jp/list/pathway")

  summed <- data.frame(data.table::fread(url, header=FALSE))
  namec <- data.frame(data.table::fread(url2, header=FALSE))
  namec <- namec$V2 %>% setNames(namec$V1)
  summed <- summed[grepl("ko", summed$V1),]
  ## No global map
  summed <- summed[!grepl("ko011", summed$V1),]
  summed <- summed[!grepl("ko012", summed$V1),]

  allpath <- unique(summed$V1)
  pathdf <- do.call(rbind, lapply(allpath, function(i) {
    tmp <- summed[summed$V1==i, ]
    int <- length(intersect(row.names(dat), tmp$V2))
    if (int>1) {
      tmpsum <- apply(dat[intersect(row.names(dat), tmp$V2),], 2, sum)
      return(c(i, tmpsum))
    } else if (int==1) {
      return(c(i, dat[intersect(row.names(dat), tmp$V2),]))
    } else {
      return(NULL)
    }
  })) %>% data.frame()
  row.names(pathdf) <- pathdf[,1]
  pathdf[,1] <- NULL
  colnames(pathdf) <- as.character(paste0("factor",seq_len(ncol(pathdf))))
  pathdf <- dplyr::mutate_all(pathdf, as.numeric)  
  if (tss) {
    pathdf <- apply(pathdf, 2, function(x) x / sum(x))
  }
  row.names(pathdf) <- namec[gsub("path:ko","map",row.names(pathdf))]
  return(pathdf)
}
