#'
#' plotCoverage
#' 
#' plot the scatters of mean coverage across groups
#' for MIDAS1
#' 
#' @param stana stana object
#' @param species candidate species
#' @param cl named list of clusters
#' @param pointSize scatter point size
#' @import ggplot2
#' 
plotCoverage <- function(stana, species, cl,
	pointSize=5) {
	if(stana@type!="MIDAS1"){stop("currently MIDAS1 only")}
	midas_merge_dir <- stana@mergeDir
	snps <- read.table(paste0(midas_merge_dir,"/",species,"/snps_summary.txt"),
                         sep="\t",header=1,row.names=1)
	ids <- row.names(snps)
	for (nm in names(cl)){
		ids[ids %in% cl[[nm]]] <- nm
	}
	snps$group <- ids
	snps$ids <- row.names(snps)
	ggplot(snps, aes(x=group, y=mean_coverage,
		fill=group)) +
	geom_jitter(shape=21, size=pointSize)+
	scale_fill_manual(values=stana@colors)+
	theme_minimal()
}

