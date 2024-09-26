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
#' @export
#' @return ggplot2 object
#' 
plotCoverage <- function(stana, species, cl=NULL,
	pointSize=5) {
	if (is.null(cl)) { cl <- stana@cl}
	midas_merge_dir <- stana@mergeDir

	if (dim(stana@snpsSummary)[1]==0) {
		if(stana@type!="MIDAS1") {
			snps <- read.table(paste0(midas_merge_dir,"/",species,"/snps_summary.txt"),
		                         sep="\t",header=1,row.names=1)			
		} else {
		    filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
    		snps <- read.table(filePath, header=1)
		}
	} else {
		snps <- stana@snpsSummary
	}
	if (stana@type=="MIDAS2") {
		snps <- snps |> dplyr::filter(species_id==species)
		row.names(snps) <- snps$sample_name
	}
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

