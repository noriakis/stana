#' checkProfile
#' 
#' Assess profile for species and return filtered species 
#' based on the number of samples for each category.
#' For MIDAS only.
#'
#' @param midas_merge_dir output directory of merge_midas.py
#' @param cl named list of sample IDs
#' @param filtNum The species with number above this threshold
#'                for each category is returned
#' @import GetoptLong
#' @export
checkProfile <- function(midas_merge_dir, cl, filtNum=2) {
  clearSn <- NULL
  clearGn <- NULL

  dirLs <- list.files(midas_merge_dir)
  specNames <- NULL
  for (d in dirLs) {
    if (dir.exists(paste0(midas_merge_dir,"/",d))){
      specNames <- c(specNames, d)
    }
  }
  freqtblSn <- NULL
  freqtblGn <- NULL
  for (sp in specNames){
    pnum <- c(sp)
    qqcat("@{sp}\n")
    qqcat("  Snps\n")
    cont <- list.files(paste0(midas_merge_dir,"/",sp))
    if ("snps_freq.txt" %in% cont) {
      snps <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_freq.txt"),
                         sep="\t",header=1,row.names=1)
      grBoolSn <- list()
      for (nm in names(cl)){
        grProfile <- length(intersect(cl[[nm]],colnames(snps)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolSn[[nm]] <- grProfile > filtNum
        pnum <- c(pnum, grProfile)
      }
      if (sum(unlist(grBoolSn))==length(cl)){
        qqcat("    @{sp} cleared filtering threshold in SNV\n")
        clearSn <- c(clearSn, sp)
      }
      freqtblSn <- rbind(freqtblSn, pnum)
    }

    pnum <- c(sp)
    if ("genes_presabs.txt" %in% cont) {
      qqcat("  Genes\n")
      genes <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_presabs.txt"),
                         sep="\t",header=1,row.names=1)
      grBoolGn <- list()
      for (nm in names(cl)){
        grProfile <- length(intersect(cl[[nm]],colnames(genes)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolGn[[nm]] <- grProfile > filtNum
        pnum <- c(pnum, grProfile)
      }
      if (sum(unlist(grBoolGn))==length(cl)){
        qqcat("    @{sp} cleared filtering threshold in genes\n")
        clearGn <- c(clearGn, sp)
      }
      freqtblGn <- rbind(freqtblGn, pnum)
    }
  }
  freqtblSn <- data.frame(freqtblSn)
  freqtblGn <- data.frame(freqtblGn)
  colnames(freqtblSn) <- c("species",names(cl))
  colnames(freqtblGn) <- c("species",names(cl))
  row.names(freqtblSn) <- seq_len(nrow(freqtblSn))
  row.names(freqtblGn) <- seq_len(nrow(freqtblGn))
  qqcat("Overall, @{length(clearSn)} species met criteria\n")
  qqcat("Overall, @{length(clearGn)} species met criteria\n")
  return(list(clearSnps=clearSn,clearGenes=clearGn,
    freqTblSn=freqtblSn, freqTblGn=freqtblGn))
}