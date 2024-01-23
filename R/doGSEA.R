#' doGSEA
#' @export
doGSEA <- function(stana, candSp=NULL, cl=NULL, eps=1e-2, how=mean,
    zeroPerc=0) {
    if (is.null(candSp)) {candSp <- stana@ids[1]}
    if (is.null(cl)) {cl <- stana@cl}
    if (length(cl)!=2) {stop("Only the two group is supported")}

    aa <- cl[[1]]
    bb <- cl[[2]]

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
    ko_df_filt <- data.frame(ko_df_filt[!rowSums(ko_df_filt==0)>dim(ko_df_filt)[2] * zeroPerc, ]) |>
        `colnames<-`(colnames(ko_df_filt))


    ## eps values
    ko_df_filt <- ko_df_filt + 1e-2
    ko_sum <- log2(apply(ko_df_filt[,intersect(inSample, aa)],1,mean) /
                   apply(ko_df_filt[,intersect(inSample, bb)],1,mean))
    ## Perform GSEA (it will take time)?
    ko_sum <- ko_sum[order(ko_sum, decreasing=TRUE)]

    ## KO to PATHWAY mapping
    url <- "https://rest.kegg.jp/link/ko/pathway"
    kopgsea <- data.frame(data.table::fread(url, header = FALSE, sep = "\t"))
    kopgsea <- kopgsea |> dplyr::filter(startsWith(V1, "path:ko"))
    kopgsea$V1 <- kopgsea$V1 |> strsplit(":") |>
      vapply("[",2,FUN.VALUE="a")
    enr <- clusterProfiler::GSEA(ko_sum,
        TERM2GENE = kopgsea, pvalueCutoff=1)
    return(enr)
}

#' calcKO
#' @export
calcKO <- function(stana, candSp=NULL, how=mean) {
    if (is.null(candSp)) {candSp <- stana@ids[1]}
    if (is.null(stana@eggNOG[[candSp]])) {stop("Please provide list of path to annotation file by `setAnnotation` function.")}
    ko_df_filt <- summariseAbundance(stana, sp = candSp,
        checkEGGNOG(annot_file=stana@eggNOG[[candSp]], "KEGG_ko"),
        how=how)
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