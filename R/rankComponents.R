#' rankComponents
#' 
#' Returns the rank of the components in the graph
#' based on the statistical values inside the graph.
#' 
#' @param stana stana object
#' @param candSp candidate species
#' @param pid KEGG pathway ID
#' @param cl if NULL, grouping in stana is used
#' @param eps pseudovalue added if log is taken
#' @param how how to combine multiple variables, default to sum
#' @param rankMethod how to rank genes
#' @param statHow how to aggregate the multiple statistical values
#' @export
rankComponents <- function(stana, pid, candSp=NULL, cl=NULL, eps=1e-2, how=sum,
    zeroPerc=0, rankMethod="modt", statHow=mean) {
    if (is.null(candSp)) {candSp <- stana@ids[1]}
    if (is.null(cl)) {cl <- stana@cl}
    if (length(cl)!=2) {stop("Only the two group is supported")}

    if (is.null(stana@kos[[candSp]])) {
        ko_tbl <- summariseAbundance(stana,sp = candSp,
                    checkEGGNOG(annot_file=stana@eggNOG[[candSp]], "KEGG_ko"),
                    how=how)
        stana@kos[[candSp]] <- ko_tbl
    } else {
        ko_tbl <- stana@kos[[candSp]]
        cat_subtle("# Using pre-computed KO table\n")
    }

    ## eps values
    cat_subtle("# ", names(cl)[1], " / ", names(cl)[2], "\n", sep="")
    ko_sum <- L2FC(ko_tbl, cl[[1]], cl[[2]], method=rankMethod, eps=eps)
    ## KO to PATHWAY mapping
    gg <- ggkegg::pathway(pid)
    gg <- gg %>% mutate(stat=ggkegg::node_numeric(ko_sum))

    types <- gg %N>% pull(type)
    ind <- which(types=="compound")


    rankVals <- lapply(ind, function(i) {
        candidate_reac <- gg %E>% dplyr::filter(to == i | from == i) %>% dplyr::pull(reaction)
        vals <- gg %N>% dplyr::filter(reaction %in% candidate_reac) %>% dplyr::pull(stat)
        do.call(statHow, list(x=vals))
    }) %>% unlist()

    gg %N>% dplyr::filter(type=="compound") %>% data.frame() %>%
        mutate(rank=rankVals) %>% arrange(desc(rank)) %>%
        select(name, rank) %>% na.omit()
}
