
#' drawEGGNOG
#' 
#' draw the graph representing the eggNOG annotation
#' 
#' @param annot_file annotation file
#' @param geneIDs geneIDs
#' @param candPlot vector consisting of 
#' "cog","ec","ko","kegg","module","reaction"
#' @export
#' 
drawEGGNOG <- function(annot_file, geneIDs, candPlot) {
  retList <- list()
  egng <- checkEGGNOG(annot_file,"all", geneIDs)
  # delete map* entry
  egng <- egng[!grepl("map",egng[,3]),]
  ap <- NULL
  for (i in geneIDs) {
    tmp <- egng[egng[,1] %in% i,]
    tmp <- tmp[ tmp[,2] %in% candPlot,]
    for (cnd in candPlot) {
      remCnd <- candPlot[ !candPlot %in% cnd ]
      if (dim(tmp[ tmp[,2] %in% cnd, ])[1]!=0) {
        for (no in tmp[ tmp[,2] %in% cnd, ][,3]) {
          for (no2 in tmp[ tmp[,2] %in% remCnd,][,3]){
            ap <- rbind(ap, c(no, no2))
          }
        }
      }
    }  
  }
  counter <- table(egng[,3])
  g <- graph_from_data_frame( ap, directed=FALSE)
  category <- NULL
  for (nn in names(V(g))) {
    category <- c(category,
      unique(as.character(egng[egng[,3]==nn,][,2])))
  }
  V(g)$category <- category
  V(g)$size <- counter[names(V(g))]
  retList[["graph"]] <- g
  plt <- ggraph(g, layout="nicely")+
            geom_edge_diagonal()+
            geom_node_point(aes(fill=.data$category, size=.data$size),shape=21)+
            geom_node_text(aes(label=.data$name, color=.data$category),
                           check_overlap=TRUE, repel=TRUE,
                           bg.color = "white", segment.color="black",
                           bg.r = .15)+
            scale_size(range=c(3,6))+
            theme_graph()
  retList[["plot"]] <- plt
  retList
}




#' checkEGGNOG
#' 
#' When the interesting genes (e.g. differentially abundant between category) are found,
#' one can query eggNOG-mapper for its functionality.
#' This function assumes one performs eggNOG-mapper to centroids.ffn of candidate species, and
#' parse the annotation results.
#' 
#' @param annot_file path to eggnog-mapper annotation
#' @param ret default to all
#' @param checkIDs only return functions related to IDs
#' @export
#' 
checkEGGNOG <- function(annot_file, ret="all", checkIDs=NULL) {
    # annot <- readr::read_table(file, comment="#", col_names=FALSE)
  annotMap <- list()
  cnmap <- data.frame(fnc=c("cog","ec","ko","kegg","module","reaction"),
    cn=c(5,11,12,13,14,15))
  con = file(annot_file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (!grepl("#", line)){
      tmp <- unlist(strsplit(line, "\t"))
      if (!is.null(checkIDs)) {
        if (tmp[1] %in% checkIDs) {
          if (ret=="all") {
            if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
            if (tmp[11]!="-"){annotMap[[tmp[1]]][["ec"]]<-tmp[11]}
            if (tmp[12]!="-"){annotMap[[tmp[1]]][["ko"]]<-tmp[12]}
            if (tmp[13]!="-"){annotMap[[tmp[1]]][["kegg"]]<-tmp[13]}
            if (tmp[14]!="-"){annotMap[[tmp[1]]][["module"]]<-tmp[14]}
            if (tmp[15]!="-"){annotMap[[tmp[1]]][["reaction"]]<-tmp[15]}
          } else if (ret=="cog") {
            if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
          } else {
            cn = as.numeric(subset(cnmap, cnmap$fnc==ret)$cn)
            if (tmp[cn]!="-"){annotMap[[tmp[1]]][[ret]]<-tmp[cn]}
          }
        } else {
          next
        }
      } else {
        if (ret=="all") {
          if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
          if (tmp[11]!="-"){annotMap[[tmp[1]]][["ec"]]<-tmp[11]}
          if (tmp[12]!="-"){annotMap[[tmp[1]]][["ko"]]<-tmp[12]}
          if (tmp[13]!="-"){annotMap[[tmp[1]]][["kegg"]]<-tmp[13]}
          if (tmp[14]!="-"){annotMap[[tmp[1]]][["module"]]<-tmp[14]}
          if (tmp[15]!="-"){annotMap[[tmp[1]]][["reaction"]]<-tmp[15]}
        } else if (ret=="cog") {
          if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
        } else {
          cn = as.numeric(subset(cnmap, cnmap$fnc==ret)$cn)
          if (tmp[cn]!="-"){annotMap[[tmp[1]]][[ret]]<-tmp[cn]}
        }
      }
    }
  }
  close(con)
  if (ret=="all"){
    # return(annotMap)
    annotDf <- NULL
    for (q in names(annotMap)){
      for (fn in c("cog","ec","ko","kegg","module","reaction")) {
        if (!is.null(annotMap[[q]][[fn]])){
          spl <- unlist(strsplit(annotMap[[q]][[fn]],","))
          for (s in spl){
            annotDf <- rbind(annotDf, c(q, fn, s))
          }
        }
      }
    }
    annotDf <- data.frame(annotDf)
    if (dim(annotDf)[1]==0) {return(NULL)}
    colnames(annotDf) <- c("id","function","functionID")
    return(annotDf)
  } else {
    annotDf <- NULL
    for (q in names(annotMap)){
      # if (!is.null(annotMap[[q]][[ret]])){
        if (ret!="cog"){
          spl <- unlist(strsplit(annotMap[[q]][[ret]],","))
        } else {
          spl <- annotMap[[q]][[ret]]
        }
        for (s in spl){
            annotDf <- rbind(annotDf, c(q, ret, s))
          }
        }
    }
    annotDf <- data.frame(annotDf)
    if (dim(annotDf)[1]==0) {return(NULL)}
    colnames(annotDf) <- c("id","function","functionID")
    return(annotDf)
  }