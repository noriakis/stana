#' checkEGGNOG
#' 
#' When the interesting genes (e.g. differentially abundant between category) are found,
#' one can query eggNOG-mapper for its functionality.
#' This function assumes one performs eggNOG-mapper to centroids.ffn of candidate species, and
#' parse the annotation results.
#' 
#' @param annot_file path to eggnog-mapper annotation
#' @param ret default to all
#' @export
#' 
checkEGGNOG <- function(annot_file, ret="all") {
    # annot <- readr::read_table(file, comment="#", col_names=FALSE)
  annotMap <- list()
  con = file(annot_file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (!grepl("#", line)){
      tmp <- unlist(strsplit(line, "\t"))
      if (grepl("COG",tmp[5])){annotMap[[tmp[1]]][["cog"]]<-tmp[5]}
      if (tmp[11]!="-"){annotMap[[tmp[1]]][["ec"]]<-tmp[11]}
      if (tmp[12]!="-"){annotMap[[tmp[1]]][["ko"]]<-tmp[12]}
      if (tmp[13]!="-"){annotMap[[tmp[1]]][["kegg"]]<-tmp[13]}
      if (tmp[14]!="-"){annotMap[[tmp[1]]][["module"]]<-tmp[14]}
      if (tmp[15]!="-"){annotMap[[tmp[1]]][["reaction"]]<-tmp[15]}
    }
  }
  close(con)
  if (ret=="all"){
    return(annotMap)
  } else {
    annotDf <- NULL
    for (q in names(annotMap)){
      if (!is.null(annotMap[[q]][[ret]])){
        spl <- unlist(strsplit(annotMap[[q]][[ret]],","))
        for (s in spl){
            annotDf <- rbind(annotDf, c(q, s))
          }
        }
    }
    annotDf <- data.frame(annotDf)
    colnames(annotDf) <- c("id","function")
    return(annotDf)
  }
}