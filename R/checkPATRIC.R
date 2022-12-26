#' checkPATRIC
#'
#' Obtain gene information from PATRIC server.
#' Input named list of genes, and returns queried results,
#' Count for functionality after removing duplicate entries,
#' and network representation of EC and KEGG pathway.
#' For midas_db_v1.2 only.
#'
#' @param genes named list of genes
#' @param whichToCount ec_description, ec_number, pathway_name, pathway_id
#' @param delSize functions with count above this threshold is included
#' @param showSubset in graph, show this number of text ordered by frequency (default: 5)
#' @param colText how to color text (default: cat, EC and KEGG)
#' @param rel use relative frequency
#' @param lyt ggraph layout (default: nicely)
#' @param nodeSize 'count' or 'degree'
#' @param skipGraph skip graph plotting
#' @import igraph ggraph BiocFileCache RCurl ggplot2
#' @export
checkPATRIC <- function(genes,
                        whichToCount="ec_description",
                        delSize=0,
                        showSubset=5,
                        colText="cat",
                        rel=FALSE,
                        lyt="nicely",
                        nodeSize="count",
                        skipGraph=FALSE) {
  allGenes <- unlist(genes)
  patricIDs <- unique(paste0(data.frame(strsplit(allGenes,"\\."))[1,],
                             ".",
                             data.frame(strsplit(allGenes,"\\."))[2,]))
  annotDf <- list()
  
  bfc <- BiocFileCache()
  qqcat("Obtaining annotations of @{length(patricIDs)} genomes\n")
  for (i in patricIDs){
    qqcat("  Obtaining information on @{i}\n")
    url <- paste0("ftp://ftp.patricbrc.org/genomes/",i,"/",i,".PATRIC.pathway.tab")
    path <- bfcrpath(bfc, url)
    tmp <- data.table::fread(path)
    tmp <- data.frame(tmp)
    annotDf[[i]] <- tmp    
  }
  annotList <- list()
  for (clname in names(genes)){
    qqcat("Checking results on cluster @{clname}\n")
    annot <- c()
    for (tmpDf in annotDf){
      greped <- tmpDf[grep(paste(genes[[clname]],collapse="|"),tmpDf$patric_id),]
      annot <- rbind(annot,
                     greped[,c("patric_id",
                               "ec_number",
                               "ec_description",
                               "pathway_id",
                               "pathway_name")]
      )
    }
    qqcat("  total of @{dim(annot)[1]} annotation obtained\n")
    annotList[[clname]][["DF"]] <- annot
    
    qqcat("  remove duplicate based on @{whichToCount}\n")
    removeDup <- annot[!duplicated(annot[,c("patric_id",
                                            whichToCount)]),]
    qqcat("  total of @{dim(removeDup)[1]} annotation obtained after removal of duplication\n")
    annotList[[clname]][["REMOVEDUP"]] <- removeDup
    collapseRes <- table(removeDup[,whichToCount])
    sorted <- collapseRes[order(collapseRes, decreasing=TRUE)]
    annotList[[clname]][["SORTED"]] <- sorted
    
    ## Make graph
    if (!skipGraph){
      ec_suff="description"
      ecget <- paste0("ec_",ec_suff)
      kegg_suff="name"
      keggget <- paste0("pathway_",kegg_suff)
      qqcat("Making graph on @{keggget} and @{ecget}\n")
      
      relec <- table(annot[!duplicated(annot[,c("patric_id",
                                       ecget)]),][,ecget])
      relkegg <- table(annot[!duplicated(annot[,c("patric_id",
                                       keggget)]),][,keggget])
      if (rel){
        relec <- relec / sum(relec)
        relkegg <- relkegg / sum(relkegg)
        counter <- c(relec, relkegg)
      } else {
        counter <- c(relec, relkegg)
      }
      qqcat("  subsetting to @{showSubset} label on each category\n")
      showNode <- c(relec[order(relec, decreasing = TRUE)][1:showSubset],
      relkegg[order(relkegg, decreasing = TRUE)][1:showSubset])
          
      # if (delOne){
      #   inc <- names(counter[counter!=1])
      #   graphInput <- annot[annot$ec_description %in% inc |
      #                         annot$pathway_name %in% inc,]
      # } else {
      #   graphInput <- annot
      # }
      
      g <- graph_from_data_frame(annot[,c(ecget,keggget)],
                                 directed = FALSE)      
      g <- simplify(g)
      
      annotList[[clname]][["GRAPH"]] <- g
      cat <- c()
      for (nm in names(V(g))){
        if (nm %in% annot[ecget][,1]){
          cat <- c(cat, "EC")
        } else {
          cat <- c(cat, "KEGG")
        }
      }
      V(g)$category <- cat
      if (nodeSize=="degree"){
        V(g)$size <- degree(g)
      } else {
        V(g)$size <- counter[names(V(g))]
      }
      
      V(g)$showText <- names(V(g)) %in% names(showNode)
      # V(g)$name <- stringr::str_wrap(V(g)$name,10)
      if (colText!="cat"){
        gp <- ggraph(g, layout=lyt)+
          geom_edge_diagonal()+
          geom_node_point(aes(size=size,fill=category),shape=21)+
          geom_node_text(aes(filter=size > delSize & showText,
                             label=name,size=size,color=size),
                         check_overlap=TRUE, repel=TRUE,
                         bg.color = "white", segment.color="black",
                         bg.r = .15)+
          scale_color_gradient(low="blue",high="red",guide="none")+
          scale_size(range=c(3,6))+
          scale_fill_manual(values=c("tomato","steelblue"),
                            name="Category")+
          theme_graph()
      } else {
        catcol <- c("tomato","steelblue")
        names(catcol) <- c("EC","KEGG")
        gp <- ggraph(g, layout=lyt)+
          geom_edge_diagonal()+
          geom_node_point(aes(size=size,fill=category),shape=21)+
          geom_node_text(aes(filter=size > delSize & showText,
                             label=name,
                             size=size,
                             color=category),
                         check_overlap=TRUE, repel=TRUE,
                         bg.color = "white", segment.color="black",
                         bg.r = .15)+
          scale_color_manual(values=catcol,
                             guide="none")+
          scale_size(range=c(3,6))+
          scale_fill_manual(values=catcol,
                            name="Category")+
          theme_graph()
        
      }
      gp <- gp + ggtitle(clname)
      annotList[[clname]][["PLOT"]] <- gp
    }
  }
  annotList
}
