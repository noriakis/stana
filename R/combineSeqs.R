#' combineSeqs
#' 
#' (experimental) use common SNV across multiple datasets to
#' infer the MSA, the same species and same database must be used.
#' 
#' @param stana_list stana list
#' @param species species ID
#' @param argList args to passed to consensusSeq
#' @param output_seq output just the sequence
#' @param minor_align minor allele alignment
#' @export
#' @return new stana object
combineSeqs <- function(stana_list, species, argList=list(), output_seq=FALSE,
    minor_align=FALSE) {
  if (!is.list(stana_list)) {stop("Please provide list of stana object")}
  each_ID <- lapply(stana_list, function(x) {
    if (minor_align) {
        ids <- paste0(row.names(x@snpsInfo[[species]]), ":minor", x@snpsInfo[[species]]$minor_allele)
    } else {
        ids <- x@snpsInfo[[species]] |> row.names()    
    }
    return(ids)
  })
  intersected <- Reduce(intersect, each_ID)
  if (length(intersected)==0) {stop("No common SNVs")}
  
  qqcat("Common SNVs: @{length(intersected)}\n")
  argList[["return_mat"]] <- TRUE
  if (minor_align) {
      intersected <- strsplit(intersected, ":minor") %>% lapply("[", 1) %>% unlist()
  }
  argList[["site_list"]] <- intersected
  
  passList <- list()
  passList[["argList"]] <- argList
  resultList <- lapply(stana_list, function(x) {
    passList[["stana"]] <- x
    passList[["species"]] <- species
    do.call(consensusSeq, passList)
  })
  
  resultList <- do.call(rbind, resultList)
  
  allele_list <- apply(resultList, 1, function(x) paste0(x, collapse=""))
  if (output_seq) {return(allele_list)}
  faName <- paste0(species,"_consensus_merged.fasta")
  if (output_seq) {
    qqcat("  Outputting consensus sequence to @{faName}\n")    	
  }
  
  fileConn<-file(faName, open="w")
  for (sample in row.names(resultList)) {
    cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
    cat(allele_list[sample], file=fileConn, append=TRUE, sep="\n")
  }
  close(fileConn)
  tre <- read.phyDat(faName, format = "fasta")
  new_stana <- new("stana")
  new_stana@fastaList[[species]] <- tre
  
  ## Metadata
  cls <- lapply(stana_list, function(x) x@cl)
  cls <- do.call(c, cls)
  
  ## Warning if same label
  ovlgr <- length(Reduce(intersect, lapply(stana_list, function(x) names(x@cl))))
  if (ovlgr>0) {qqcat("Duplicate label found in group\n")}
  new_stana@ids <- species
  new_stana@type <- stana_list[[1]]@type
  new_stana@cl <- cls
  new_stana@colors <- getColors(cls)
  if (!output_seq) {
    unlink(faName)
  }
  return(new_stana)
}

