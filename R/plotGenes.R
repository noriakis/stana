#'
#' plotGenes
#' 
#' plot the violin plot for specific genes
#' 
#' @param stana stana object
#' @param species candidate species
#' @param geneID gene ID to be plotted
#' @param cl named list of clusters
#' @import ggplot2
#' @export
#' @return ggplot
#' 
plotGenes <- function(stana, species, geneID, cl=NULL) {
	if (is.null(cl)) { cl <- stana@cl}
	geneDf <- stana@genes[[species]]
	if (intersect(geneID,row.names(geneDf)) |> length() > 0) {
		geneDf <- geneDf[intersect(geneID,row.names(geneDf)), ]
	} else {
		stop("geneID not in data.frame")
	}
	geneDf$geneID <- row.names(geneDf)
	df <- geneDf |> tidyr::pivot_longer(1:ncol(geneDf)-1)
	ids <- df$name
	for (nm in names(cl)){
		ids[ids %in% cl[[nm]]] <- nm
	}
	df$group <- ids
	ggplot(df, aes(x=group, y=value,
		fill=group)) +
	geom_violin()+
	geom_jitter(shape=21, size=2)+
	scale_fill_manual(values=stana@colors)+
	facet_wrap(.~geneID)+
	theme_minimal()
}
