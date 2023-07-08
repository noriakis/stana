#'
#' compareGenes
#' 
#' compare the gene abundances by recursively performing Wilcoxon rank-sum tests
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
compareGenes <- function(stana, species, geneID=NULL, cl=NULL, argList=list(),
	verbose_zero=FALSE) {
	if (is.null(cl)) { cl <- stana@cl}
	midas_merge_dir <- stana@mergeDir
	geneDf <- stana@genes[[species]]
	incSamples <- intersect(colnames(geneDf), unlist(cl))
	geneDf <- geneDf[, incSamples]
	if (!is.null(geneID)) {
		if (intersect(geneID,row.names(geneDf)) |> length() > 0) {
			geneDf <- geneDf[intersect(geneID,row.names(geneDf)), ]
		} else {
			stop("geneID not in data.frame")
		}		
	}
	if (length(names(cl))!=2) {stop("Only two grouping is allowed")}
	resList <- list()
	qqcat("Testing total of @{nrow(geneDf)}\n")
	for (i in row.names(geneDf)) {
		x <- geneDf[i, intersect(colnames(geneDf), cl[[1]])] |> as.numeric()
		y <- geneDf[i, intersect(colnames(geneDf), cl[[2]])] |> as.numeric()
		if (verbose_zero) {
			if (sum(x)==0 ) {qqcat("  In @{i}, Group @{names(cl)[1]} all zero\n")}
			if (sum(y)==0 ) {qqcat("  In @{i}, Group @{names(cl)[2]} all zero\n")}			
		}
		argList[["x"]] <- x
		argList[["y"]] <- y
		argList[["paired"]] <- FALSE
		resList[[i]] <- do.call(exactRankTests::wilcox.exact, argList)
	}
	resList
}
