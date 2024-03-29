#'
#' plotDist
#' 
#' Perform distance calculation based on the matrix
#' stored in stana, and plot the heatmap with grouping
#' information
#' 
#' @param stana stana object
#' @param sp species to be examined
#' @param cl named list of samples
#' @param target tree, snps, genes, kos, fasta
#' @param AAfunc if choose `fasta`, provide function for calculating distance
#' @param AAargs provided to `AAfunc`
#' @param distMethod distance method passed to dist() (default, manhattan)
#' @param distArg passed to `dist()
#' @export
plotDist <- function(stana, sp, cl=NULL, AAfunc=dist.ml, AAargs=list(),
    target="snps", distMethod="manhattan", distArg=list()) {
      if (is.null(cl)) {cl <- stana@cl}
      if (length(sp)>1) {stop("Please specify one species")}
        cat_subtle("# Performing dist in ", sp, " target is ", target, "\n", sep="")
        if (target=="snps") {
            snps <- stana@snps[[sp]]
            if (!is.null(stana@includeSNVID[[sp]])) {
                cat_subtle("# The set SNV ID information (", length(stana@includeSNVID[[sp]]), ") is used.\n")
                snps <- snps[stana@includeSNVID[[sp]], ]
            }
            ## Replace zero depth with NA
            snps[snps == -1] <- NA
        } else if (target=="tree") {
            if (!is.null(stana@treeList[[sp]])) {
              tre <- stana@treeList[[sp]]
              ## cophenetic distance
              d <- as.dist(ape::cophenetic.phylo(tre))
              sn <- attr(d, "Labels")
            } else {
              stop("No tree found in stana@treeList")
            }
        } else if (target=="fasta") {
          AAargs[["x"]] <- stana@fastaList[[sp]]
          d <- do.call(AAfunc, AAargs)
          sn <- attr(d, "Labels")
        } else if (target=="kos") {
            snps <- stana@kos[[sp]]          
        } else {
            snps <- stana@genes[[sp]]
        }

        if (!(target %in% c("tree","fasta"))) {
          distArg[["x"]] <- t(snps)
          distArg[["method"]] <- distMethod
          d <- do.call(dist, distArg)
          sn <- attr(d, "Labels")      
        }
        met <- data.frame(listToNV(cl))
        met <- data.frame(met[sn, ]) %>% `row.names<-`(sn)
        pheatmap(d, annotation_row=met)
}