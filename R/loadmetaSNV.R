


#' loadmetaSNV
#' 
#' Assess and store profile for species for metaSNV.
#'
#' @param metasnv_out_dir output directory of merge_midas.py
#' @param cl named list of sample IDs
#' @param just_species just return species id
#' @param candSp candidate species ID
#' @import GetoptLong
#' @export
loadmetaSNV <- function(metasnv_out_dir, cl=NULL,
                        just_species=FALSE, candSp=NULL) {
  stana <- new("stana")
  stana@type <- "metaSNV"
  snpList <- list()
  stana@mergeDir <- metasnv_out_dir
  dirLs <- list.files(metasnv_out_dir)
  if ("filtered" %in% dirLs) {
    if (dir.exists(paste0(metasnv_out_dir,"/filtered/pop"))){
      freqList <- list.files(paste0(metasnv_out_dir,"/filtered/pop"))
      spList <- unlist(lapply(strsplit(freqList, ".filtered"),"[",1))
      if (just_species) {return(spList)}
      stana@ids <- spList
      if (!is.null(candSp)) {spList <- candSp}
      for (sp in spList) {
        qqcat("  Loading @{sp}\n")
        df <- read.table(paste0(metasnv_out_dir,"/filtered/pop/",sp,".filtered.freq"),
          header=1, row.names=1)
        snpList[[sp]] <- df
      }
    }
  }
  if (!is.null(cl)) {stana@cl <- cl}
  stana <- initializeStana(stana,cl)
  stana@snps <- snpList
  stana
}
