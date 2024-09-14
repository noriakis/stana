
#' calcGF
#' @param stana stana object
#' @param candSp candidate species ID
#' @param how how to summarize multiple gene CN assigned to the same KO
#' @param annot eggNOG or manual
#' @param column When eggNOG, which family to summarize, default to KEGG_ko
#' @export
calcGF <- function(stana, candSp=NULL, how=sum, annot="eggNOG", column="KEGG_ko") {
	checkID(stana, candSp)
    if (is.null(candSp)) {cat("Species not specified, the first ID will be used:", stana@ids[1]);
        candSp <- stana@ids[1]
    }
    if (annot=="eggNOG") {
	    if (is.null(stana@eggNOG[[candSp]])) {stop("Please provide list of path to annotation file by `setAnnotation` function.")}
        cat_subtle("# Using eggNOG slot\n")
	    ko_df_filt <- summariseAbundance(stana, sp = candSp,
	        checkEGGNOG(annot_file=stana@eggNOG[[candSp]], column),
	        how=how)  	
    } else {
    	if (is.null(stana@map[[candSp]])) {stop("Please set mapping data.frame in map slot using `setMap`.")}
        cat_subtle("# Using map slot\n")
	    ko_df_filt <- summariseAbundance(stana, sp = candSp,
	        stana@map[[candSp]],
	        how=how)
    }
    stana@kos[[candSp]] <- ko_df_filt
    stana
}