

#' plotHeatmap
#' @description Plot a heatmap of gene copy number.
#' @details Plot a heatmap with functional annotation by simplifyEnrichment.
#' The annotations are provided by `eggNOG` slot or `map` slot. If loaded type is MIDAS,
#' The function automatically fetches the functional annotation from PATRIC server.
#' @param stana stana object
#' @param sp candidate species
#' @param cl grouping named list
#' @param k k-means param for splitting gene
#' @param geneID gene id to subset
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
#' @param filter_max_value max value for copy numbers (default to 50), coupled with filter_max_frac
#' @param variable If specified other than 0, subset to top-{variable} variation gene numbers.
#' @importFrom ComplexHeatmap Heatmap
#' @export
#' 
plotHeatmap <- function(stana, sp, cl=NULL, k=10, mat=NULL, seed=1,
	geneID=NULL, variable=0,
	fnc="KEGG_Pathway", removeHigh=TRUE, removeAdditional=NULL, max_words=10,
    filter_zero_frac=0.8, filter_max_frac=0, filter_max_value=50) {
	set.seed(seed)

	if (!is.null(mat)) {
		df <- mat
	} else {
		df <- stana@genes[[sp]]
        ## Filter
    	if (!is.null(geneID)) {
            cat_subtle("# Ignoring filtering options\n")
			df <- df[intersect(row.names(df), geneID), ]
		} else {
            df <- df[!rowSums(df == 0) > ncol(df) * filter_zero_frac,]
            df <- df[!rowSums(df > filter_max_value) > ncol(df) * filter_max_frac,]                     
        }
	}


    cat_subtle("# In resulting matrix, max: ", max(df), " min: ", min(df), "\n", sep="")
	cat_subtle("# Dimension: ", dim(df)[1], " ", dim(df)[2], "\n", sep="")
	
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
    
    if (variable!=0) {
        tmpvar <- matrixStats::rowVars(expr %>% as.matrix()) %>% sort(decreasing=TRUE) %>%
            head(variable) %>% names()
        expr <- expr[tmpvar, ]
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
			      argList=list(max_words = max_words, fontsize_range=c(10,20)))
			)
		return(hm)		
	} else {
		cat_subtle("# Looking for the annotation file by eggNOG-mapper v2\n")
		if (is.null(stana@eggNOG[[sp]])) {stop("No annotation file provided to slot `eggNOG`")}
		cat_subtle("# Loading annotation\n")
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
			      argList=list(max_words = max_words, fontsize_range=c(10,20)))
			)
		return(hm)
	}

}