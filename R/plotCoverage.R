#'
#' plotCoverage
#' 
#' plot the scatters of mean coverage across groups
#' 
#' @param midas_merge_dir output of midas merge
#' @param species candidate species
#' @param cl named list of clusters
#' @import ggplot2
#' 
plotCoverage <- function(midas_merge_dir, species, cl,
	pointSize=5) {
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
	theme_minimal()
}

