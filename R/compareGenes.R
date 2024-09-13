#'
#' @title compareGenes
#' 
#' @description Compare the gene abundances by recursively performing Wilcoxon rank-sum tests
#' Use exactRankTests package, and multiple comparison adjustment is not performed
#' 
#' @param stana stana object
#' @param species candidate species
#' @param geneID gene ID to be plotted
#' @param cl named list of clusters
#' @param argList argument list for wilcox.exact
#' @param verbose_zero show zero abundance genes
#' @import ggplot2
#' @export
#' @return ggplot
#' 
compareGenes <- function(stana,
	species=NULL, geneID=NULL,
	cl=NULL, argList=list(),
	verbose_zero=FALSE) {
	
	if (is.null(cl)) { cl <- stana@cl}
	if (length(cl)>2) {stop("Two group is supported")}
	if (is.null(species)) {species <- stana@clearGenes[1]}
	midas_merge_dir <- stana@mergeDir
	geneDf <- stana@genes[[species]]
	incSamples <- intersect(colnames(geneDf), unlist(cl))
	geneDf <- geneDf[, c("gene_id", ..incSamples)]
	if (!is.null(geneID)) {
		if (intersect(geneID, geneDf$gene_id) %>% length() > 0) {
			geneDf <- geneDf[ gene_id %in% intersect(geneID, geneDf$gene_id), ]
		} else {
			stop("geneID not in data.frame")
		}		
	}
	if (length(names(cl))!=2) {stop("Only two grouping is allowed")}
	resList <- list()
	cat_subtle("# Testing total of ", nrow(geneDf), "\n", sep="")
	for (i in geneDf$gene_id) {
		x <- geneDf[gene_id == i, intersect(colnames(geneDf), cl[[1]]), with=FALSE] %>% unlist() %>% as.numeric()
		y <- geneDf[gene_id == i, intersect(colnames(geneDf), cl[[2]]), with=FALSE] %>% unlist() %>% as.numeric()
		if (verbose_zero) {
			if (sum(x)==0 ) {cat_subtle("  In ",i, "Group ",names(cl)[1]," all zero\n")}
			if (sum(y)==0 ) {cat_subtle("  In ",i, "Group ",names(cl)[2]," all zero\n")}			
		}
		argList[["x"]] <- x
		argList[["y"]] <- y
		argList[["paired"]] <- FALSE
		resList[[i]] <- do.call(exactRankTests::wilcox.exact, argList)
	}
	resList
}
