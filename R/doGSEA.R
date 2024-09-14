#' doGSEA
#' 
#' Based on KEGG database, GSEA will be performed by the function.
#' clusterProfiler::GSEA is used.
#' By default this uses moderated t-statistics for ranking of the genes.
#' 
#' @param stana stana object
#' @param candSp candidate species
#' @param cl if NULL, grouping in stana is used
#' @param eps pseudovalue added if log is taken
#' @param how how to combine multiple variables, default to sum
#' @param rankMethod how to rank genes
#' @param target pathway or module
#' @param zeroPerc genes >= the percentage of count zero sample will be excluded.
#' Default to zero, not recommended in GSEA
#' @param bg_filter filter the background for those in table
#' @param background background TERM2GENE data
#' @param gseaArgs list of argument passed to GSEA or fgsea function
#' @importFrom MKmisc mod.t.test
#' @return GSEA results from clusterProfiler
#' @export
doGSEA <- function(stana, candSp=NULL, cl=NULL, eps=1e-2, how=sum,
    zeroPerc=0, rankMethod="modt", target="pathway", bg_filter=TRUE,
    background=NULL, gseaArgs=list()) {
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
        stana@kos[[candSp]] <- ko_df_filt
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
    if (is.null(gseaArgs[["TERM2GENE"]])) {
        if (!is.null(background)) {
            gseaArgs[["TERM2GENE"]] <- background
        } else {
            gseaArgs[["TERM2GENE"]] <- kopgsea
        }        
    }
    if (is.null(gseaArgs[["pvalueCutoff"]])) {
        ## Return all the value regardless of P
        gseaArgs[["pvalueCutoff"]] <- 1    
    }
    
    gseaArgs[["geneList"]] <- ko_sum
    enr <- do.call(clusterProfiler::GSEA, gseaArgs)

    if (!is.null(stana@gsea[[candSp]])) {
        cat_subtle("# Overriding previous GSEA result\n")
    }
    stana@gsea[[candSp]] <- enr
    return(stana)
}

#' plotGSEA
#' 
#' plot GSEA results
#' 
#' @param stana stana object or list of stana object
#' @param dataset_names dataset names
#' @param padjThreshold threshold of adjusted p-values
#' @param return_graph return only the tbl_graph
#' @param layout graph layout in ggraph
#' @export
#' @return ggplot object
plotGSEA <- function(stana, dataset_names=NULL, padjThreshold=0.05,
    return_graph=FALSE, layout="kk") {
    if (is.list(stana)) {
        num <- length(stana)
    } else {
        num <- 1
        stana <- list(stana)
    }
    if (is.null(dataset_names)) {
        dataset_names <- paste0("dataset_",seq_len(num))
    }

    dat <- do.call(rbind, lapply(seq_along(stana), function(e) {
        gsea_res <- stana[[e]]@gsea
        summed <- do.call(rbind, lapply(names(gsea_res), function(x) {
            tmp <- subset(gsea_res[[x]]@result, p.adjust<padjThreshold)
            if (dim(tmp)[1]!=0) {
                tmp[["dataset"]] <- dataset_names[e]
                tmp[["species_name"]] <- x
                return(tmp)
            } else {
                return(NULL)
            }
        }))
    }))

    if (is.null(dat)) {
        stop("There is no significant results")
    }

    raw <- dat %>%
        as_tibble() %>%
        dplyr::select(ID, NES, p.adjust, species_name, dataset)

    raw2 <- bind_rows(raw %>% select(ID, species_name, dataset) %>% `colnames<-`(c("from","to","dataset")))
    raw2 <- na.omit(raw2)
    gg <- as_tbl_graph(raw2)

    if (return_graph) {
        return(gg)
    }

    gg <- gg %N>%
        mutate(pathway=startsWith(name, "ko")) %>%
        mutate(species=!pathway) %>%
        mutate(degree=centrality_degree(mode="all"))


    ggraph(gg, layout=layout) +
      geom_edge_parallel(aes(color=dataset)) +
      geom_node_point(aes(filter=pathway, size=degree), shape=15) +
      geom_node_text(aes(label=name, filter=pathway, size=degree), repel=TRUE, bg.colour="white") +
      geom_node_point(aes(filter=species, size=degree)) +
      geom_node_text(aes(label=name, filter=species, size=degree), repel=TRUE, bg.colour="white") +
      theme_graph()+
      scale_size(range=c(2,6))

}


#' addGeneAbundance
#' 
#' Add the specified gene copy number to metadata
#' 
#' @param stana stana object
#' @param candSp candidate species
#' @param IDs gene IDs to add
#' @param target KO or genes
#' @param how how to combine multiple IDs
#' @param newCol new column name
#' @param discNumeric convert discrete value to numeric
#' @param disc discretize the abundance by the threshold. function for calculating
#' threshold, like {median}
#' @param convert conversion such as log10
#' @export
#' @return stana
addGeneAbundance <- function(stana, candSp, IDs,
    target="KO", how=sum, newCol="gene",
    disc=NULL, discNumeric=TRUE, convert=NULL) {
    if (target=="KO") {
        ints <- intersect(row.names(stana@kos[[candSp]]), IDs)
        subMat <- stana@kos[[candSp]][ints, ]
    } else {
        ints <- intersect(row.names(stana@genes[[candSp]]), IDs)
        subMat <- stana@genes[[candSp]][ints, ]
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
    if (!is.null(convert)) {
        adda <- do.call(convert, list(x=adda))
    }
    meta[[newCol]] <- adda[row.names(stana@meta)]
    stana@meta <- meta
    return(stana)
}

#' reverse_annot
#' @param stana stana object
#' @param candSp candidate species
#' @param candidate candidate gene family ids
#' @param col column to use in eggNOG mapper annotation
#' @export
reverseAnnot <- function(stana, candSp, candidate, col="KEGG_ko") {
    annot <- checkEGGNOG(annot_file=stana@eggNOG[[candSp]], col)
    annot %>% dplyr::filter(.data[["value"]] %in% candidate) %>%
        dplyr::pull(ID) %>% unique()
}


#' obtain_positions
#' @param stana stana object
#' @param candSp candidate species
#' @param geneID gene IDs
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