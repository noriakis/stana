
#' drawEGGNOG
#' 
#' draw the graph representing the eggNOG annotation
#' 
#' @param annot_file annotation file
#' @param geneIDs geneIDs
#' @param candPlot vector consisting of 
#' "COG_category" "Description" "Preferred_name" "GOs"           
#' "EC"             "KEGG_ko"        "KEGG_Pathway"   "KEGG_Module"    "KEGG_Reaction" 
#' "KEGG_rclass"    "BRITE"          "KEGG_TC"        "CAZy"           "BiGG_Reaction" 
#' "PFAMs"
#' @export
#' 
drawEGGNOG <- function(annot_file, geneIDs, candPlot) {
  retList <- list()
  egng <- checkEGGNOG(annot_file,"all", geneIDs)
  # delete map* entry (in KEGG)
  egng <- egng[!grepl("map",egng$value),]
  print(egng)
  ap <- NULL
  for (i in geneIDs) {
    tmp <- egng[egng$ID == i,]
    tmp <- tmp[ tmp$name %in% candPlot,]
    for (cnd in candPlot) {
      remCnd <- candPlot[ !candPlot %in% cnd ]
      if (dim(tmp[ tmp$name %in% cnd, ])[1]!=0) {
        for (no in tmp[ tmp$name %in% cnd, ]$value) {
          for (no2 in tmp[ tmp$name %in% remCnd,]$value){
            ap <- rbind(ap, c(no, no2))
          }
        }
      }
    }  
  }

  counter <- table(egng$value)
  g <- graph_from_data_frame( ap, directed=FALSE)
  category <- NULL
  for (nn in names(V(g))) {
    category <- c(category,
      unique(as.character(egng[egng$value==nn,]$name)))
  }
  V(g)$category <- category
  V(g)$size <- as.numeric(counter[names(V(g))])
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
#' 
#' Choices of IDs in `ret`:
#' "COG_category" "Description" "Preferred_name" "GOs"           
#' "EC"             "KEGG_ko"        "KEGG_Pathway"   "KEGG_Module"    "KEGG_Reaction" 
#' "KEGG_rclass"    "BRITE"          "KEGG_TC"        "CAZy"           "BiGG_Reaction" 
#' "PFAMs"
#' 
#' @param annot_file path to eggnog-mapper annotation
#' @param ret default to all
#' @param checkIDs only return functions related to IDs
#' @param fill data.table::fread argument
#' @importFrom data.table fread
#' @export
#' 
checkEGGNOG <- function(annot_file, ret="all", checkIDs=NULL, fill=TRUE) {
  ann <- data.table::fread(annot_file,
                           skip=4,sep="\t",fill=fill)
  if (ret!="all") {
    parsed <- ann[,c("#query",ret), with=FALSE] |>
      tidyr::pivot_longer(-1) |>
      dplyr::filter(value!="-") |>
      dplyr::mutate(value = strsplit(as.character(value), ",")) |>
      tidyr::unnest(value)
  } else {
    parsed <- ann[, !c("evalue","score"), with=FALSE] |>
      tidyr::pivot_longer(-1) |>
      dplyr::filter(value!="-") |>
      dplyr::mutate(value = strsplit(as.character(value), ",")) |>
      tidyr::unnest(value)
  }
  parsed <- parsed |> `colnames<-`(c("ID","name","value"))
  if (!is.null(checkIDs)) {
    parsed <- parsed |> 
      dplyr::filter(ID %in% checkIDs)
  }
  parsed






  ## Below is an old function
  #   # annot <- readr::read_table(file, comment="#", col_names=FALSE)
  # annotMap <- list()
  # cnmap <- data.frame(fnc=c("cog","ec","ko","kegg","module","reaction"),
  #   cn=c(5,11,12,13,14,15))
  # con = file(annot_file, "r")
  # while ( TRUE ) {
  #   line = readLines(con, n = 1)
  #   if ( length(line) == 0 ) {
  #     break
  #   }
  #   if (!grepl("#", line)){
  #     tmp <- unlist(strsplit(line, "\t"))
  #     if (!is.null(checkIDs)) {
  #       if (tmp[1] %in% checkIDs) {
  #         if (ret=="all") {
  #           if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
  #           if (tmp[11]!="-"){annotMap[[tmp[1]]][["ec"]]<-tmp[11]}
  #           if (tmp[12]!="-"){annotMap[[tmp[1]]][["ko"]]<-tmp[12]}
  #           if (tmp[13]!="-"){annotMap[[tmp[1]]][["kegg"]]<-tmp[13]}
  #           if (tmp[14]!="-"){annotMap[[tmp[1]]][["module"]]<-tmp[14]}
  #           if (tmp[15]!="-"){annotMap[[tmp[1]]][["reaction"]]<-tmp[15]}
  #         } else if (ret=="cog") {
  #           if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
  #         } else {
  #           cn = as.numeric(subset(cnmap, cnmap$fnc==ret)$cn)
  #           if (tmp[cn]!="-"){annotMap[[tmp[1]]][[ret]]<-tmp[cn]}
  #         }
  #       } else {
  #         next
  #       }
  #     } else {
  #       if (ret=="all") {
  #         if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
  #         if (tmp[11]!="-"){annotMap[[tmp[1]]][["ec"]]<-tmp[11]}
  #         if (tmp[12]!="-"){annotMap[[tmp[1]]][["ko"]]<-tmp[12]}
  #         if (tmp[13]!="-"){annotMap[[tmp[1]]][["kegg"]]<-tmp[13]}
  #         if (tmp[14]!="-"){annotMap[[tmp[1]]][["module"]]<-tmp[14]}
  #         if (tmp[15]!="-"){annotMap[[tmp[1]]][["reaction"]]<-tmp[15]}
  #       } else if (ret=="cog") {
  #         if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
  #       } else {
  #         cn = as.numeric(subset(cnmap, cnmap$fnc==ret)$cn)
  #         if (tmp[cn]!="-"){annotMap[[tmp[1]]][[ret]]<-tmp[cn]}
  #       }
  #     }
  #   }
  # }
  # close(con)
  # if (ret=="all"){
  #   # return(annotMap)
  #   annotDf <- NULL
  #   for (q in names(annotMap)){
  #     for (fn in c("cog","ec","ko","kegg","module","reaction")) {
  #       if (!is.null(annotMap[[q]][[fn]])){
  #         spl <- unlist(strsplit(annotMap[[q]][[fn]],","))
  #         for (s in spl){
  #           annotDf <- rbind(annotDf, c(q, fn, s))
  #         }
  #       }
  #     }
  #   }
  #   annotDf <- data.frame(annotDf)
  #   if (dim(annotDf)[1]==0) {return(NULL)}
  #   colnames(annotDf) <- c("id","function","functionID")
  #   return(annotDf)
  # } else {
  #   annotDf <- NULL
  #   for (q in names(annotMap)){
  #     # if (!is.null(annotMap[[q]][[ret]])){
  #       if (ret!="cog"){
  #         spl <- unlist(strsplit(annotMap[[q]][[ret]],","))
  #       } else {
  #         spl <- annotMap[[q]][[ret]]
  #       }
  #       for (s in spl){
  #           annotDf <- rbind(annotDf, c(q, ret, s))
  #         }
  #       }
  #   }
  #   annotDf <- data.frame(annotDf)
  #   if (dim(annotDf)[1]==0) {return(NULL)}
  #   colnames(annotDf) <- c("id","function","functionID")
  #   return(annotDf)
  }

  
#' summariseAbundance
#' 
#' Given tibble with `ID` and `value` column,
#' Obtain data.frame of genes consisting of `ID` within `value`
#' and summarise them by default `mean`.
#' 
#' @param stana stana object
#' @param sp candidate species
#' @param anno annotation tibble obtained by checkEGGNOG
#' @param how summarising function, default to mean
#' @return data.frame
#' @export
#' 
summariseAbundance <- function(stana, sp, anno, how=mean) {
  geneDf <- stana@genes[[sp]]
  merged <- list()
  annoval <- anno$value |> unique()
  for (i in annoval) {
    candID <- (anno |> dplyr::filter(anno$value==i))$ID
    ints <- intersect(row.names(geneDf), candID)
    if (length(ints)>0) {
      merged[[i]] <- apply(geneDf[ints,], 2, how)
    }
  }
  do.call(rbind, merged)
}