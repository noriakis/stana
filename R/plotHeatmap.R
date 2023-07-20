

#' plotHeatmap
#' 
#' plot a heatmap with functional annotation by simplifyEnrichment
#' Typically, MIDAS and MIDAS2 output are used
#' @param stana stana object
#' @param sp candidate species
#' @param cl grouping named list
#' @param km k-means param for splitting gene
#' @param mat if customized matrix is to be used. Row will be gene names
#' and column sample names
#' @param fnc One of `KEGG_Pathway` or `KEGG_Module`, when eggNOG annotation is used
#' @importFrom ComplexHeatmap Heatmap
#' @export
#' 
plotHeatmap <- function(stana, sp, cl=NULL, km=10, mat=NULL, seed=1,
	fnc="KEGG_Pathway") {
	set.seed(seed)

	if (!is.null(mat)) {
		df <- mat
	} else {
		df <- stana@genes[[sp]]
	}
	df <- df |> head(500)
	qqcat("Dimension: @{dim(df)[1]}, @{dim(df)[2]}\n")
	if (is.null(cl)) {cl <- stana@cl}
	if (length(cl)==0) {cl <- list("all"=colnames(df))}

	expr <- df[,intersect(colnames(df),as.character(unlist(cl)))]

	ord <- NULL
	spl <- NULL
	for (nm in names(cl)){
	    inc <- intersect(colnames(expr), cl[[nm]])
	    ord <- c(ord, inc)
	    spl <- c(spl, rep(nm, length(inc)))
	}

	km = kmeans(expr, centers = 20)$cluster

	if (stana@type=="MIDAS") {
	    hm <- Heatmap(expr, row_split = km,
	        column_split=spl,
	        show_column_names = FALSE,border=TRUE,name="Copy number",
	        show_row_names = FALSE, show_row_dend = FALSE,
	        col=c("steelblue","white","tomato")) + 
			rowAnnotation(
			    keywords = stana::anno_PATRIC_keywords(split = km,
			      genes = rownames(expr), fnc="pathway_name",
			      removeHigh=TRUE, 
			      argList=list(max_words = 10))
			)
		return(hm)		
	} else if (stana@type=="MIDAS2") {
		qqcat("MIDAS2, looking for the annotation file by eggNOG-mapper v2\n")
		if (is.null(stana@eggNOG[[sp]])) {stop("No annotation file provided to slot `eggNOG`")}
		qqcat("Loading annotation\n")
		tib <- checkEGGNOG(stana@eggNOG[[sp]], ret=fnc)
	    if (fnc=="KEGG_Pathway") {
     		mp <- data.table::fread("https://rest.kegg.jp/list/pathway", header=FALSE)
     		mp$value <- paste0("ko",mp$V1 |> strsplit("map") |> sapply("[",2))
	    } else if (fnc=="KEGG_Module") {
     		mp <- data.table::fread("https://rest.kegg.jp/list/pathway", header=FALSE)
     		mp$value <- mp$V1
	    } else {return(1)}
	    tib <- merge(tib, mp, by="value")
	    hm <- Heatmap(expr, row_split = km,
	        column_split=spl,
	        show_column_names = FALSE,border=TRUE,name="Copy number",
	        show_row_names = FALSE, show_row_dend = FALSE,
	        col=c("steelblue","white","tomato")) + 
			rowAnnotation(
			    keywords = stana::anno_eggNOG_keywords(split = km,
			      genes = rownames(expr), tib=tib,
			      removeHigh=TRUE, 
			      argList=list(max_words = 10))
			)
		return(hm)
	} else {}


}