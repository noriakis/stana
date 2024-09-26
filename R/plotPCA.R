#'
#' plotPCA
#' 
#' plot the PCA results from snp or gene table
#' 
#' @param stana output of midas merge
#' @param target snps, genes
#' @param species candidate species
#' @param cl named list of clusters
#' @param pointSize scatter point size
#' @param ltype line type for ellipse
#' @param useBlend use ggblend for point
#' @param ignoreZero ignore the SNVs containing -1
#' @param replaceZero replace the SNVs -1 (zero depth) with zero
#' @import ggplot2
#' @importFrom stats prcomp
#' @export
#' @return ggplot
#' 
plotPCA <- function(stana, species, cl=NULL, target="snps",
	ignoreZero=FALSE, replaceZero=FALSE,
	pointSize=5, ltype=2, useBlend=FALSE) {
	cols <- stana@colors
	pcaList <- list()
	if (is.null(cl)) {cl <- stana@cl}

	for (sp in species) {
		
		if (target=="snps") {
			df <- stana@snps[[sp]]
			if (!is.null(stana@includeSNVID[[sp]])) {
				cat_subtle("# The set SNV ID information (", length(stana@includeSNVID[[sp]]), ") is used.\n")
				df <- df[stana@includeSNVID[[sp]], ]
			}
		} else if (target=="genes") {
			df <- stana@genes[[sp]]
		} else if (target=="kos") {
			df <- stana@kos[[sp]]	
		} else {
			stop("please specify snps, genes or kos")
		}
		
    	if (length(cl)==0) {} else {
			df <- df[,intersect(colnames(df),
				as.character(unlist(cl)))]    		
    	}
		if (ignoreZero) {
			df <- df[rowSums(df==-1)==0,]
		}
		if (replaceZero) {
			## [CAUTION]
			df[df==-1] <- 0
		}
		qqcat("After filtering: @{dim(df)[1]} SNVs\n")
		if (length(cl)!=0) {
			ids <- colnames(df)
			for (nm in names(cl)){
				ids[ids %in% cl[[nm]]] <- nm
			}
		} else {
			ids <- "No_Group"
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
		    theme_minimal()+ggtitle(sp)
		} else {
		    p <- pcdf |>
		    ggplot(aes(x=PC1, y=PC2,
		    	color=Group, fill=Group))+
		    geom_point(size=pointSize,
		    	alpha=0.7)+
			scale_color_manual(values=cols)+
			scale_fill_manual(values=cols)+
			stat_ellipse(lty=ltype)+
		    theme_minimal()+ggtitle(sp)
		}
		pcaList[[sp]] <- p
	}
	pcaList
}

