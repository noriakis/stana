

#' plotHeatmap
#' 
#' plot a heatmap with functional annotation by simplifyEnrichment
#' Typically, MIDAS and MIDAS2 output are used
#' @param stana stana object
#' @param sp candidate species
#' @param cl grouping named list
#' @param k k-means param for splitting gene
#' @param mat if customized matrix is to be used. Row will be gene names
#' and column sample names
#' @param fnc One of `KEGG_Pathway` or `KEGG_Module`, when eggNOG annotation is used
#' @param removeHigh remove high frequent words (preset)
#' @param removeAdditional remove additional words specified
#' @param max_words max words to plot
#' @param seed random seed
#' @param filter_zero_frac genes with zero abundance over fraction of samples as this value
#' are removed before sample filtering. As typically gene matrix is large, for further filtering, please use `mat` option
#' @param filter_max_frac remove genes with values below `filter_max_value` in this fraction of sample
#' @param filter_max_value max value for copy numbers
#' @importFrom ComplexHeatmap Heatmap
#' @export
#' 
plotHeatmap <- function(stana, sp, cl=NULL, k=10, mat=NULL, seed=1,
	fnc="KEGG_Pathway", removeHigh=TRUE, removeAdditional=NULL, max_words=10,
    filter_zero_frac=0.8, filter_max_frac=Inf, filter_max_value=5) {
	set.seed(seed)

	if (!is.null(mat)) {
		df <- mat
	} else {
		df <- stana@genes[[sp]]
        ## Filter
        df <- df[!rowSums(df == 0) > ncol(df) * filter_zero_frac,]
        df <- df[!rowSums(df > filter_max_value) > ncol(df) * filter_max_frac,]
	}

    qqcat("In resulting matrix, max: @{max(df)}, min: @{min(df)}\n")

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

	km = kmeans(expr, centers = k)$cluster

	if (stana@type=="MIDAS1") {
	    hm <- Heatmap(expr, row_split = km,
	        column_split=spl,
	        show_column_names = FALSE,border=TRUE,name="Copy number",
	        show_row_names = FALSE, show_row_dend = FALSE,
	        col=c("steelblue","white","tomato")) + 
			rowAnnotation(
			    keywords = stana::anno_PATRIC_keywords(split = km,
			      genes = rownames(expr), fnc="pathway_name",
			      removeHigh=removeHigh,  removeAdditional=removeAdditional,
			      argList=list(max_words = max_words))
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
     		mp <- data.table::fread("https://rest.kegg.jp/list/module", header=FALSE)
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
			      removeHigh=removeHigh,  removeAdditional=removeAdditional,
			      argList=list(max_words = max_words))
			)
		return(hm)
	} else {
		stop("Currently MIDAS2 and MIDAS are supported.")
	}


}