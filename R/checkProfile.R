
#' loadInStrain
#' 
#' For InStrain, you should provide `compare` output
#' with --bams options implmented from version 1.6.
#' It computes MAF of biallelic position.
#' Filters of popSNV or conSNV will be added.
#' 
#' @param compare_out_dir output directory of compare.
#' need `output` in the directory
#' @param sample_threshold SNVs with at least this number of samples are included
#' @param candidate_species candidate species ID, e.g. `GUT_GENOME000024`
#' @param just_species if TRUE, returns the vector of species IDs
#' 
#' @export
loadInStrain <- function(compare_out_dir, sample_threshold, candidate_species, just_species=FALSE) {
    output_list <- list.files(paste0(compare_out_dir,"/output"))
    stana <- new("stana")
    stana@type <- "InStrain"
    stana@mergeDir <- compare_out_dir
    snps <- list()
    if (sum(grepl("pooled_SNV_data",output_list))!=0) {
        path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data.tsv.gz",output_list)])
        tbl <- data.table::fread(path)
        keys_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data_keys.tsv",output_list)])
        keys <- data.table::fread(keys_path)
        info_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_info.tsv.gz",output_list)])
        info <- data.table::fread(info_path)
        sps <- unique(paste0(sapply(strsplit(keys$scaffold, "_"),"[",1),"_",sapply(strsplit(keys$scaffold, "_"),"[",2)))
        if (just_species) {
            return(sps)
        }
        stana@ids <- sps
        qqcat("Candidate species: @{candidate_species}\n")
        candKeys <- keys[grepl(candidate_species,keys$scaffold),]$key
        qqcat("  Candidate key numbers: @{length(candKeys)}\n")
        ret <- function(x) {
            smp <- keys[keys$key==x[1],"sample"]$sample
            scaff <- keys[keys$key==x[2],"scaffold"]$scaffold
            det <- c(smp, paste0(scaff,"|",x[3]))
            majmin <- x[4:7][x[4:7]!=0][1:2]
            if (length(unique(majmin)!=1)) {
                major <- names(majmin)[which.max(majmin)]
                minor <- names(majmin)[which.min(majmin)]
                det <- c(det, major, minor,
                         as.numeric((x[4:7] / sum(x[4:7]))[minor]))
            } else {
                ## Alphabetical order
                m <- sort(names(majmin))
                det <- c(det, m[1], m[2],
                         as.numeric((x[4:7] / sum(x[4:7]))[m[2]]))
            }
            det
        }
        mafs <- NULL

        sdt <- tbl[tbl$scaffold %in% candKeys,]
        qqcat("  Dimension of pooled SNV table for species: @{dim(sdt)[1]}")
        # Only bi-allelic position
        sdt <- sdt[apply(sdt[,4:7],1,function(x)sum(x==0))==2,]
        maf <- apply(sdt, 1, function(x) ret(x))
        maf <- data.frame(t(maf))

        maf <- maf |> `colnames<-`(c("id","scaffold_position","major_allele","minor_allele","maf"))
        maf$scaffold_position <- paste0(maf$scaffold_position,"|",maf$major_allele,">",maf$minor_allele)
        sample_snv <- tidyr::pivot_wider(maf, id_cols=id, names_from = scaffold_position, values_from = maf)
        thresh <- apply(sample_snv[2:ncol(sample_snv)], 2, function(x) sum(!is.na(x))) > sample_threshold
        snps[[candidate_species]] <- sample_snv[,c("id",names(thresh[thresh]))]
        stana@snps <- snps
    }
    return(stana)
}


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
loadmetaSNV <- function(metasnv_out_dir, cl=NULL,
                        filtType="group",
                        filtNum=2, filtPer=0.8) {
  stana <- new("stana")
  stana@type <- "metaSNV"
  snpList <- list()
  stana@mergeDir <- metasnv_out_dir
  dirLs <- list.files(metasnv_out_dir)
  if ("filtered" %in% dirLs) {
    if (dir.exists(paste0(metasnv_out_dir,"/filtered/pop"))){
      freqList <- list.files(paste0(metasnv_out_dir,"/filtered/pop"))
      spList <- unlist(lapply(strsplit(freqList, ".filtered"),"[",1))
      stana@ids <- spList
      for (sp in spList) {
        df <- read.table(paste0(metasnv_out_dir,"/filtered/pop/",sp,".filtered.freq"),
          header=1, row.names=1)
        snpList[[sp]] <- df
      }
    }
  }
  stana@snps <- snpList
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
#' @param candSp candidate species ID
#' @param geneType "presabs" or "copynum"
#' @import GetoptLong
#' @export
loadMIDAS <- function(midas_merge_dir,
  cl, filtType="group", candSp=NULL,
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

  if (filtType=="whole") {
    checkCl <- list("whole"=as.character(unlist(cl)))
  } else if (filtType=="group") {
    checkCl <- cl
  } else {
    stop("Please specify whole or group")
  }
  stana@sampleFilter=filtType

  if (!is.null(candSp)) {specNames <- candSp}
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
      for (nm in names(checkCl)){
        grProfile <- length(intersect(checkCl[[nm]],colnames(snps)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolSn[[nm]] <- grProfile > filtNum | grProfile > length(checkCl[[nm]])*filtPer
        pnum <- c(pnum, grProfile)
      }
      if (sum(unlist(grBoolSn))==length(checkCl)){
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
      for (nm in names(checkCl)){
        grProfile <- length(intersect(checkCl[[nm]],colnames(genes)))
        qqcat("    @{nm} @{grProfile}\n")
        grBoolGn[[nm]] <- grProfile > filtNum | grProfile > length(checkCl[[nm]])*filtPer
        pnum <- c(pnum, grProfile)
      }
      if (sum(unlist(grBoolGn))==length(checkCl)){
        qqcat("    @{sp} cleared filtering threshold in genes\n")
        clearGn <- c(clearGn, sp)
      }
      freqtblGn <- rbind(freqtblGn, pnum)
    }
  }
  freqtblSn <- data.frame(freqtblSn)
  freqtblGn <- data.frame(freqtblGn)
  colnames(freqtblSn) <- c("species",names(checkCl))
  colnames(freqtblGn) <- c("species",names(checkCl))
  row.names(freqtblSn) <- seq_len(nrow(freqtblSn))
  row.names(freqtblGn) <- seq_len(nrow(freqtblGn))
  stana@clearSnps <- clearSn
  stana@clearGenes <- clearGn
  stana@freqTableSnps <- freqtblSn
  stana@freqTableGenes <- freqtblGn
  stana@snps <- snpList
  stana@genes <- geneList
  
  stana <- initializeStana(stana,cl)

  qqcat("Overall, @{length(clearSn)} species met criteria in SNPs\n")
  qqcat("Overall, @{length(clearGn)} species met criteria in genes\n")
  stana
}

#' getColors
#' @import RColorBrewer
#' @noRd
getColors <- function(cl){
  numgr <- length(names(cl))
  cols <- brewer.pal(numgr, "PuOr") 
}

#' initializeStana
#' @noRd
initializeStana <- function(stana,cl) {
  stana@colors <- getColors(cl)
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
#' @param candSp candidate species ID
#' @param taxtbl tax table, row.names: 6-digits MIDAS2 ID and `GTDB species` column.
#' @export
#' 
loadMIDAS2 <- function(midas_merge_dir,
                        cl,
                        filtNum=2,
                        filtPer=0.8,
                        taxtbl=NULL,
                        candSp=NULL,
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

  if (filtType=="whole") {
    checkCl <- list("whole"=as.character(unlist(cl)))
  } else if (filtType=="group") {
    checkCl <- cl
  } else {
    stop("Please specify whole or group")
  }
  stana@sampleFilter=filtType

  qqcat("SNPS\v")
  if (!is.null(candSp)) {specNames <- candSp} else {
    specNames <- list.files(paste0(midas_merge_dir,"/snps"))
  }
  for (i in specNames){
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
          for (clName in names(checkCl)){
              int <- intersect(colnames(df), checkCl[[clName]])
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(cl[[clName]])*filtPer)
          }
          if (sum(checkPass)==length(names(checkCl))){
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
  if (!is.null(candSp)) {specNames <- candSp} else {
    specNames <- list.files(paste0(midas_merge_dir,"/genes"))
  }
  for (i in specNames){
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
          for (clName in names(checkCl)){
              int <- intersect(colnames(df), checkCl[[clName]])
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(cl[[clName]])*filtPer)
          }
          if (sum(checkPass)==length(names(checkCl))){
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
  stana <- initializeStana(stana,cl)
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
#' For MIDAS only. Deprecated, use load* functions instead.
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
