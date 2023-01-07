#'
#' plotPCA
#' 
#' plot the PCA results from snp or gene table
#' 
#' @param staba output of midas merge
#' @param species candidate species
#' @param cl named list of clusters
#' @param pointSize scatter point size
#' @param cols colors to point
#' @param ltype line type for ellipse
#' @param useBlend use ggblend for point
#' @import ggplot2
#' @import ggblend
#' @export
#' 
plotPCA <- function(stana, species, cl, target="snps",
	pointSize=5, cols=palette(), ltype=2, useBlend=FALSE) {
	pcaList <- list()
	for (sp in species) {
		if (target=="snps") {
			df <- stana@snps[[sp]]
		} else {
			df <- stana@genes[[sp]]
		}
		if (!is.null(cl)) {
			ids <- colnames(df)
			for (nm in names(cl)){
				ids[ids %in% cl[[nm]]] <- nm
			}
		}
	    pc <- prcomp(t(df))
	    pcdf <- data.frame(cbind(pc$x[,1:2], ids),
	    	check.names=FALSE)
	    pcdf[,1] <- as.numeric(pcdf[,1])
	    pcdf[,2] <- as.numeric(pcdf[,2])
	    pcdf$Group <- pcdf[,3]


	    if (useBlend){
		    p <- pcdf |>
		    ggplot(aes(x=PC1, y=PC2,
		    	color=Group, fill=Group))+
		    geom_point(size=pointSize)|> ggblend::blend("darken")+
			scale_color_manual(values=cols)+
			scale_fill_manual(values=cols)+
			stat_ellipse(lty=ltype)+
		    theme_minimal()
		} else {
		    p <- pcdf |>
		    ggplot(aes(x=PC1, y=PC2,
		    	color=Group, fill=Group))+
		    geom_point(size=pointSize)+
			scale_color_manual(values=cols)+
			scale_fill_manual(values=cols)+
			stat_ellipse(lty=ltype)+
		    theme_minimal()
		}
		pcaList[[sp]] <- p
	}
	pcaList
}

