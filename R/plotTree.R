#' plotTree
#' 
#' infer the tree and plot
#' 
#' @param stana
#' @param species
#' 
#' @export
plotTree <- function(stana, species, cl=NULL) {
	if (is.null(cl)) {cl <- stana@cl}
	tre <- stana@fastaList[[species]]
	dm <- dist.ml(tre, "F81")
	tre <- NJ(dm)
	stana@treeList[[species]] <- tre
    tre <- groupOTU(tre, cl)
    tp <- ggtree(tre,
        layout='circular',branch.length = "none") + # Return cladogram by default
        geom_tippoint(size=3, aes(color=.data$group)) +
        ggtitle(species)+
        scale_color_manual(values=stana@colors)
    stana@treePlotList[[species]] <- tp
    stana
}