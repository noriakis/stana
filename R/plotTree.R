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
#' @param tree_args passed to dist function in phangorn
#' @param dist_method dist method in phangorn, default to dist.ml
#' @param branch.length branch length, default to "none", cladogram
#' 
#' @export
plotTree <- function(stana, species=NULL, cl=NULL,
	dist_method="dist.ml",
	tree_args=list(), branch.length="none") {
	if (is.null(cl)) {cl <- stana@cl}
	if (is.null(species)) {species <- names(stana@fastaList)}
	for (sp in species) {
		tre <- stana@fastaList[[sp]]
		tree_args[["x"]] <- tre
		dm <- do.call(dist_method, tree_args)
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