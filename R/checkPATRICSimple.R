#' checkPATRICSimple
#'
#' Obtain gene information from PATRIC server.
#' Input named list of genes, and returns named list of functions.
#' For midas_db_v1.2 only.
#'
#' @param genes named list of genes
#' @param whichToCount ec_description, ec_number, pathway_name, pathway_id
#' @import BiocFileCache RCurl
#' @export
checkPATRICSimple <- function(genes,
                        whichToCount="ec_description") {
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
    qqcat("  remove duplicate based on @{whichToCount}\n")
    removeDup <- annot[!duplicated(annot[,c("patric_id",
                                            whichToCount)]),]
    qqcat("  total of @{dim(removeDup)[1]} annotation obtained after removal of duplication\n")
    annotList[[clname]] <- as.character(removeDup[,whichToCount])
  }
  annotList
}
