#'
#' plotMAF
#' 
#' plot the box plot of MAF for specificed SNV position
#' 
#' @param stana stana object
#' @param species candidate species
#' @param SNV snv position name
#' @param cl named list of clusters
#' @param deleteZeroDepth delete position with zero depth
#' @import ggplot2
#' @export
#' @return ggplot
#' 
plotMAF <- function(stana, species, SNV, cl=NULL, deleteZeroDepth=TRUE) {
	if (is.null(cl)) { cl <- stana@cl}
	snvDf <- stana@snps[[species]]
	if (intersect(SNV, row.names(snvDf)) |> length() > 0) {
		snvDf <- snvDf[intersect(SNV,row.names(snvDf)), ]
	} else {
		stop("snv ID not in data.frame")
	}
	snvDf$snvID <- row.names(snvDf)
	df <- snvDf |> tidyr::pivot_longer(1:ncol(snvDf)-1)
	# ids <- df$name
	# for (nm in names(cl)){
	# 	ids[ids %in% cl[[nm]]] <- nm
	# }
	df$group <- listToNV(stana@cl)[df$name]
	if (deleteZeroDepth) {
		df <- df |> dplyr::filter(df$value!=-1)
	} else {
		ch <- df$value; ch[ch == -1] <- NA
		df$value <- ch
	}
	ggplot(df, aes(x=group, y=value,
		fill=group)) +
	geom_violin()+
	geom_jitter(shape=21, size=2)+
	scale_fill_manual(values=stana@colors)+
	facet_wrap(.~snvID)+
	theme_minimal()
}
