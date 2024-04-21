#'
#' plotGenes
#' 
#' @description Plot the plot for specific genes
#' @details Plot the gene copy number between groups based on the specified geom.
#' The default is to use gene copy number, but gene family abundance can also be specified.
#' 
#' @param stana stana object
#' @param species candidate species
#' @param geneID gene ID to be plotted
#' @param cl named list of clusters
#' @param target genes or KOs (gene families)
#' @param return_df return data only
#' @param geom default to geom_violin, can be changed to the geom like `geom_boxplot`
#' @import ggplot2
#' @export
#' @return ggplot
#' 
plotGenes <- function(stana, species, geneID, target="genes", cl=NULL, return_df=FALSE,
    geom=geom_violin()) {
	if (is.null(cl)) { cl <- stana@cl}
	if (target=="genes") {
    	geneDf <- stana@genes[[species]]
	} else {
		geneDf <- data.frame(stana@kos[[species]])
	}
	if (intersect(geneID,row.names(geneDf)) |> length() > 0) {
		geneDf <- geneDf[intersect(geneID,row.names(geneDf)), ]
	} else {
		stop("geneID not in data.frame")
	}
	geneDf$geneID <- row.names(geneDf)
	df <- geneDf |> tidyr::pivot_longer(1:ncol(geneDf)-1)
	df$group <- listToNV(stana@cl)[df$name]
    
    if (return_df){
        return(df)
    }

	ggplot(df, aes(x=group, y=value,
		fill=group)) +
	geom+
	geom_jitter(shape=21, size=2)+
	scale_fill_manual(values=stana@colors)+
	facet_wrap(.~geneID)+
	cowplot::theme_cowplot()+
    cowplot::panel_border()
}
