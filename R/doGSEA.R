#' doGSEA
#' 
#' Based on KEGG database, GSEA will be performed by the function.
#' clusterProfiler::GSEA is used.
#' By default this uses moderated t-statistics for ranking of the genes.
#' 
#' @param zeroPerc genes >= the percentage of count zero sample will be excluded.
#' Default to zero, not recommended in GSEA
#' @param bg_filter filter the background for those in table
#' @importFrom MKmisc mod.t.test
#' @return GSEA results from clusterProfiler
#' @export
doGSEA <- function(stana, candSp=NULL, cl=NULL, eps=1e-2, how=mean,
    zeroPerc=0, rankMethod="modt", target="pathway", bg_filter=TRUE) {
    if (is.null(candSp)) {candSp <- stana@ids[1]}
    if (is.null(cl)) {cl <- stana@cl}
    if (length(cl)!=2) {stop("Only the two group is supported")}

    aa <- cl[[1]]
    bb <- cl[[2]]
    cat(names(cl)[1], "/", names(cl)[2], "\n");
    geneMat <- stana@genes[[candSp]]
    inSample <- colnames(geneMat)
    if (length(intersect(inSample, aa))==0) {stop("No sample available")}
    if (length(intersect(inSample, bb))==0) {stop("No sample available")}

    if (is.null(stana@kos[[candSp]])) {
        cat("KO abundance not found, calculating based on annotation ...\n")
        ko_df_filt <- summariseAbundance(stana, sp = candSp,
            checkEGGNOG(annot_file=stana@eggNOG[[candSp]], "KEGG_ko"),
            how=how)
    } else {
        ko_df_filt <- stana@kos[[candSp]]
    }

    ## If filter based on number of zero per KOs
    sub <- unique(row.names(ko_df_filt))
    ko_df_filt <- data.frame(ko_df_filt[ rowSums(ko_df_filt!=0) >= (dim(ko_df_filt)[2] * zeroPerc), ]) |>
        `colnames<-`(colnames(ko_df_filt))

    ## eps values
    ko_sum <- L2FC(ko_df_filt, aa, bb, method=rankMethod, eps=eps)
    ## Perform GSEA (it will take time)?
    # print(ko_sum)
    ko_sum <- ko_sum[order(ko_sum, decreasing=TRUE)]

    ## KO to PATHWAY mapping
    
    url <- "https://rest.kegg.jp/link/ko/pathway"
    if (target!="pathway") {
	    url <- "https://rest.kegg.jp/link/ko/module"	
    }
    
    bfc <- BiocFileCache()
    path <- bfcrpath(bfc, url)
  
    kopgsea <- data.frame(data.table::fread(path, header = FALSE, sep = "\t"))
    if (target=="pathway") {
	    kopgsea <- kopgsea |> dplyr::filter(startsWith(V1, "path:ko"))	
    }
    kopgsea$V1 <- kopgsea$V1 |> strsplit(":") |>
      vapply("[",2,FUN.VALUE="a")

    if (bg_filter) {
        kopgsea <- kopgsea[kopgsea$V2 %in% sub, ]
    }
    print(dim(kopgsea))
    ## Return all the value regardless of P
    enr <- GSEA(ko_sum,
        TERM2GENE = kopgsea, pvalueCutoff=1)
    return(enr)
}


#' addGeneAbundance
#' 
#' Add the specified gene copy number to metadata
#' 
#' @param disc discretize the abundance by the threshold. function for calculating
#' threshold, like {median}
#' @export
addGeneAbundance <- function(stana, candSp, IDs,
    target="KO", how=sum, newCol="gene",
    disc=NULL, discNumeric=TRUE) {
    if (target=="KO") {
        subMat <- stana@kos[[candSp]][IDs, ]
    } else {
        subMat <- stana@genes[[candSp]][IDs, ]
    }
    if (length(IDs)>1) {
        adda <- apply(subMat, 2, how)    
    } else {
        adda <- subMat
    }
    nm <- names(adda)
    if (!is.null(disc)) {
        thresh <- do.call(disc, list("x"=adda))
        if (discNumeric) {
            adda <- as.numeric(adda > thresh)        
        } else {
            adda <- adda > thresh
        }
        names(adda) <- nm
    }
    meta <- stana@meta
    meta[[newCol]] <- adda[row.names(stana@meta)]
    stana@meta <- meta
    return(stana)
}

#' L2FC
#' 
#' Report various statistics for use in GSEA or visualization
#' 
#' @param mat row corresponds to gene, and column samples
#' @param l1 level1
#' @param l2 level2
#' @param method gmean, amean, or t
#' @param eps pseudocount added when calculating log
#' @noRd
L2FC <- function(mat, l1, l2, method="t", eps=0) {
    if (method == "gmean") {
        l1_mean <- apply(log2(mat[, intersect(colnames(mat), l1)] + eps), 1, mean)
        l2_mean <- apply(log2(mat[, intersect(colnames(mat), l2)] + eps), 1, mean)
        return(l1_mean - l2_mean)
    } else if (method == "t") {
        res <- lapply(row.names(mat), function(m) {
            tres <- t.test(mat[m, intersect(colnames(mat), l1)],
                mat[m, intersect(colnames(mat), l2)])
            as.numeric(tres$statistic)
        }) %>% unlist()
        names(res) <- row.names(mat)
        return(res)
    } else if (method == "amean" ) {
        l1_mean <- apply(mat[, intersect(colnames(mat), l1)], 1, mean)
        l2_mean <- apply(mat[, intersect(colnames(mat), l2)], 1, mean)
        return(log2((l1_mean+eps) / (l2_mean+eps)))
    } else {
    	## Moderated t.test (limma)
    	ordered.mat <- mat[, c( intersect(colnames(mat), l1), intersect(colnames(mat), l2) )]
    	gr <- c( rep("l1", length(intersect(colnames(mat), l1))), rep("l2", length(intersect(colnames(mat), l2))) )
    	res <- MKmisc::mod.t.test(as.matrix(ordered.mat), gr)
    	modt <- res$t
    	names(modt) <- row.names(res)
    	return(modt)
    }
}


#' calcKO
#' @param stana stana object
#' @param candSp candidate species ID
#' @param how how to summarize multiple gene CN assigned to the same KO
#' @param annot eggNOG or manual
#' @export
calcKO <- function(stana, candSp=NULL, how=sum, annot="eggNOG") {
	checkID(stana, candSp)
    if (is.null(candSp)) {cat("Species not specified, the first ID will be used:", stana@ids[1]);
        candSp <- stana@ids[1]
    }
    if (annot=="eggNOG") {
	    if (is.null(stana@eggNOG[[candSp]])) {stop("Please provide list of path to annotation file by `setAnnotation` function.")}
	    ko_df_filt <- summariseAbundance(stana, sp = candSp,
	        checkEGGNOG(annot_file=stana@eggNOG[[candSp]], "KEGG_ko"),
	        how=how)    	
    } else {
    	if (is.null(stana@map[[candSp]])) {stop("Please set mapping data.frame in map slot using `setMap`.")}
	    ko_df_filt <- summariseAbundance(stana, sp = candSp,
	        stana@map[[candSp]],
	        how=how)
    }

    stana@kos[[candSp]] <- ko_df_filt
    stana
}

#' reverse_annot
#' @export
reverseAnnot <- function(stana, candSp, candidate, col="KEGG_ko") {
    annot <- checkEGGNOG(annot_file=stana@eggNOG[[candSp]], col)
    annot %>% filter(.data[["value"]] %in% candidate) %>%
        pull(ID) %>% unique()
}


#' obtain_positions
#' @export
obtainPositions <- function(stana, candSp, geneID) {
    tmp <- stana@snpsInfo[[candSp]]
    tmp[tmp$gene_id %in% geneID, ] %>% row.names()
}



#' @noRd
checkID <- function(stana, candSp) {
	if (!(candSp %in% stana@ids)) {
		stop("No specified species available in this stana object.")
	}
}