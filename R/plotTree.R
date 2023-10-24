#' plotTree
#' 
#' infer the tree and plot
#' use dist.ml and NJ for inference
#' 
#' @param stana stana object
#' @param species species to plot
#' If NULL, first species in fasta list is assigned
#' @param cl optional, cluster to plot
#' @param model dist.ml model
#' @param branch.length branch length, default to "none", cladogram
#' 
#' @export
plotTree <- function(stana, species=NULL, cl=NULL,
	model="F81", branch.length="none") {
	if (is.null(cl)) {cl <- stana@cl}
	if (is.null(species)) {species <- names(stana@fastaList)}
	
	for (sp in species) {
		tre <- stana@fastaList[[sp]]
		dm <- dist.ml(tre, model)
		tre <- NJ(dm)
		stana@treeList[[sp]] <- tre
	    tre <- groupOTU(tre, cl)
	    tp <- ggtree(tre,
	        layout='circular',branch.length = branch.length) + # Return cladogram by default
	        geom_tippoint(size=3, aes(color=.data$group)) +
	        ggtitle(sp)+
	        scale_color_manual(values=stana@colors)
	    stana@treePlotList[[sp]] <- tp		
	}
    stana
}