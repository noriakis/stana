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
    colnames(annotDf) <- c("id","function","functionID")
    return(annotDf)
  }