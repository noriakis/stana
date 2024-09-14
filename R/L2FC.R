#' L2FC
#' 
#' Report various statistics for use in GSEA or visualization
#' 
#' @param mat row corresponds to gene, and column samples
#' @param l1 level1
#' @param l2 level2
#' @param method gmean, amean, or t
#' @param eps pseudocount added when calculating log
#' @return named vector of statistical values
#' @noRd
L2FC <- function(mat, l1, l2, method="t", eps=0) {
    if (method == "gmean") {
        l1_mean <- apply(log2(mat[, intersect(colnames(mat), l1)] + eps), 1, mean)
        l2_mean <- apply(log2(mat[, intersect(colnames(mat), l2)] + eps), 1, mean)
        return(l1_mean - l2_mean)
    } else if (method == "t") {
        res <- lapply(row.names(mat), function(m) {
            tres <- t.test(mat[m, intersect(colnames(mat), l1)],
                mat[m, intersect(colnames(mat), l2)])
            as.numeric(tres$statistic)
        }) %>% unlist()
        names(res) <- row.names(mat)
        return(res)
    } else if (method == "amean" ) {
        l1_mean <- apply(mat[, intersect(colnames(mat), l1)], 1, mean)
        l2_mean <- apply(mat[, intersect(colnames(mat), l2)], 1, mean)
        return(log2((l1_mean+eps) / (l2_mean+eps)))
    } else {
    	## Moderated t.test (limma)
    	ordered.mat <- mat[, c( intersect(colnames(mat), l1), intersect(colnames(mat), l2) )]
    	gr <- c( rep("l1", length(intersect(colnames(mat), l1))), rep("l2", length(intersect(colnames(mat), l2))) )
    	res <- MKmisc::mod.t.test(as.matrix(ordered.mat), gr)
    	modt <- res$t
    	names(modt) <- row.names(res)
    	return(modt)
    }
}