#' loadmetaSNV
#' 
#' Assess and store profile for species and return filtered species 
#' based on the number of samples for each category or whole population.
#' For metaSNV only.
#'
#' @param metasnv_out_dir output directory of merge_midas.py
#' @param cl named list of sample IDs
#' @import GetoptLong
#' @export
loadmetaSNV <- function(metasnv_out_dir, cl=NULL) {
  stana <- new("stana")
  stana@type <- "metaSNV"
  stana@mergeDir <- metasnv_out_dir
  dirLs <- list.files(metasnv_out_dir)
  if ("filtered" %in% dirLs) {
    if (dir.exists(paste0(metasnv_out_dir,"/filtered/pop"))){
      freqList <- list.files(paste0(metasnv_out_dir,"/filtered/pop"))
      spList <- unlist(lapply(strsplit(freqList, ".filtered"),"[",1))
      stana@ids <- spList
      for (freq in freqList) {
        df <- read.table(paste0(metasnv_out_dir,"/filtered/pop/",freq,
          header=1, row.names=1))
      }
    }
  }
  stana
}


#' loadMIDAS
#' 
#' Assess and store profile for species and return filtered species 
#' based on the number of samples for each category or whole population.
#' For MIDAS only.
#'
#' @param midas_merge_dir output directory of merge_midas.py
#' @param cl named list of sample IDs
#' @param filtType "whole" or "group"
#' @param filtNum The species with number above this threshold
#'                for each category is returned
#' @param filtPer filter by percentage
#' @param geneType "presabs" or "copynum"
#' @import GetoptLong
#' @export
loadMIDAS <- function(midas_merge_dir,
  cl, filtType="group",
  filtNum=2, filtPer=0.8,
  geneType="copynum") {
  stana <- new("stana")
  stana@type <- "MIDAS1"
  stana@mergeDir <- midas_merge_dir
  stana@sampleFilter <- filtType
  stana@sampleFilterVal <- filtNum
  stana@sampleFilterPer <- filtPer
  dirLs <- list.files(midas_merge_dir)
  specNames <- NULL
  for (d in dirLs) {
    if (dir.exists(paste0(midas_merge_dir,"/",d))){
      specNames <- c(specNames, d)
    }
  }
  stana@ids <- specNames
  stana@geneType <- geneType
  clearSn <- NULL
  clearGn <- NULL
  freqtblSn <- NULL
  freqtblGn <- NULL
  snpList <- list()
  geneList <- list()
  for (sp in specNames){
    pnum <- c(sp)
    qqcat("@{sp}\n")
    qqcat("  Snps\n")
    cont <- list.files(paste0(midas_merge_dir,"/",sp))
    if ("snps_freq.txt" %in% cont) {
      snps <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_freq.txt"),
                         sep="\t",header=1,row.names=1)
      snpList[[sp]] <- snps
      grBoolSn <- list()
      for (nm in names(cl)){
        grProfile <- length(intersect(cl[[nm]],colnames(snps)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolSn[[nm]] <- grProfile > filtNum | grProfile > length(cl[[nm]])*filtPer
        pnum <- c(pnum, grProfile)
      }
      if (sum(unlist(grBoolSn))==length(cl)){
        qqcat("    @{sp} cleared filtering threshold in SNV\n")
        clearSn <- c(clearSn, sp)
      }
      freqtblSn <- rbind(freqtblSn, pnum)
    }

    pnum <- c(sp)
    if (paste0("genes_",geneType,".txt") %in% cont) {
      qqcat("  Genes\n")
      genes <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_",geneType,".txt"),
                         sep="\t",header=1,row.names=1)
      geneList[[sp]] <- genes
      grBoolGn <- list()
      for (nm in names(cl)){
        grProfile <- length(intersect(cl[[nm]],colnames(genes)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolGn[[nm]] <- grProfile > filtNum | grProfile > length(cl[[nm]])*filtPer
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
  stana@clearSnps <- clearSn
  stana@clearGenes <- clearGn
  stana@freqTableSnps <- freqtblSn
  stana@freqTableGenes <- freqtblGn
  stana@snps <- snpList
  stana@genes <- geneList
  
  stana <- initializeStana(stana)

  qqcat("Overall, @{length(clearSn)} species met criteria in SNPs\n")
  qqcat("Overall, @{length(clearGn)} species met criteria in genes\n")
  stana
}

#' initializeStana
#' @noRd
initializeStana <- function(stana) {
  faList <- vector("list", length(stana@ids))
  names(faList) <- stana@ids
  stana@fastaList <- faList
  treeList <- vector("list", length(stana@ids))
  names(treeList) <- stana@ids
  stana@treeList <- treeList
  treePlotList <- vector("list", length(stana@ids))
  names(treePlotList) <- stana@ids
  stana@treePlotList <- treePlotList
  adonisList <- vector("list", length(stana@ids))
  names(adonisList) <- stana@ids
  stana@adonisList <- adonisList
  stana
}

#' loadMIDAS2
#' 
#' load the MIDAS2 merge output.
#' The location to lz4 binary must be added to PATH.
#' Tax table can be loaded from downloaded MIDAS2 db directory (metadata.tsv).
#' If provided, additionally show tax names.
#' @param midas_merge_dir path to merged directory
#' @param cl named list of category for samples
#' @param filtNum the species with the samples above this number will be returned
#' @param filtPer filter by percentage
#' @param taxtbl tax table, row.names: 6-digits MIDAS2 ID and `GTDB species` column.
#' @export
#' 
loadMIDAS2 <- function(midas_merge_dir,
                        cl,
                        filtNum=2,
                        filtPer=0.8,
                        taxtbl=NULL,
                        filtType="group",
                        geneType="copynum") {
  stana <- new("stana")
  stana@type <- "MIDAS2"
  stana@mergeDir <- midas_merge_dir
  stana@geneType <- geneType
  stana@sampleFilter <- filtType
  stana@sampleFilterVal <- filtNum
  stana@sampleFilterPer <- filtPer
  snpList <- list()
  geneList <- list()
  clearGn <- NULL
  clearGnSp <- NULL
  clearSn <- NULL
  clearSnSp <- NULL

  qqcat("SNPS\v")
  for (i in list.files(paste0(midas_merge_dir,"/snps"))){
      if (!is.na(as.numeric(i))) {
          qqcat("  @{i}\n")
          if (!is.null(taxtbl)){
              spnm <- taxtbl[i,]$`GTDB species`
              qqcat("  @{spnm}\n")
          }
          cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_freqs.tsv.lz4")
          cnd <- gsub(".lz4","",cnc)
          system2("lz4", args=c("-d","-f",
                                paste0(getwd(),"/",cnc),
                                paste0(getwd(),"/",cnd)),
                  stdout=FALSE, stderr=FALSE)
          df <- read.table(cnd, row.names=1, header=1)
          snpList[[i]] <- df
          qqcat("    Number of snps: @{dim(df)[1]}\n")
          qqcat("    Number of samples: @{dim(df)[2]}\n")
          # snpList[[i]] <- df
          
          checkPass <- NULL
          for (clName in names(cl)){
              int <- intersect(colnames(df), cl[[clName]])
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(cl[[clName]])*filtPer)
          }
          if (sum(checkPass)==length(names(cl))){
              qqcat("      Passed the filter\n")
              clearSn <- c(clearSn, i)
              if (!is.null(taxtbl)){
                  clearSnSp <- c(clearSnSp, spnm)
              }
          } else {
              qqcat("      Not passed the filter\n")            
          }
          
          unlink(paste0(getwd(),"/",cnd))
      }
  }
  qqcat("Genes\n")
  for (i in list.files(paste0(midas_merge_dir,"/genes"))){
      if (!is.na(as.numeric(i))) {
          qqcat("  @{i}\n")
          if (!is.null(taxtbl)){
              spnm <- taxtbl[i,]$`GTDB species`
              qqcat("  @{spnm}\n")
          }
          cnc <- paste0(midas_merge_dir,
                        "/genes/",i,"/",i,".genes_",geneType,".tsv.lz4")
          cnd <- gsub(".lz4","",cnc)
          system2("lz4", args=c("-d", "-f",
                                paste0(getwd(),"/",cnc),
                                paste0(getwd(),"/",cnd)),
                  stdout=FALSE, stderr=FALSE)
          df <- read.table(cnd, row.names=1, header=1)
          qqcat("    Number of genes: @{dim(df)[1]}\n")
          qqcat("    Number of samples: @{dim(df)[2]}\n")
          geneList[[i]] <- df
          
          checkPass <- NULL
          for (clName in names(cl)){
              int <- intersect(colnames(df), cl[[clName]])
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(cl[[clName]])*filtPer)
          }
          if (sum(checkPass)==length(names(cl))){
              qqcat("      Passed the filter\n")
              clearGn <- c(clearGn, i)
              if (!is.null(taxtbl)){
                  clearGnSp <- c(clearGnSp, spnm)
              }
          } else {
              qqcat("      Not passed the filter\n")            
          }
          unlink(paste0(getwd(),"/",cnd))
      }
  }
  stana@clearSnps <- clearSn
  stana@clearGenes <- clearGn
  stana@snps <- snpList
  stana@genes <- geneList
  stana@ids <- union(names(geneList),names(snpList))
  stana <- initializeStana(stana)
  if (!is.null(taxtbl)){
    stana@clearSnpsSpecies <- clearSnSp
    stana@clearGenesSpecies <- clearGnSp
  }
  stana
}

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
