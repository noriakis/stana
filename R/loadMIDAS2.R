#' loadMIDAS2
#' 
#' load the MIDAS2 merge command output.
#' The location to lz4 binary must be added to PATH.
#' Taxonomy table can be loaded from downloaded MIDAS2 db directory (metadata.tsv).
#' If provided, additionally show tax names.
#' 
#' @param midas_merge_dir path to merged directory
#' @param cl named list of category for samples
#' @param filtNum the species with the samples above this number will be returned
#' @param filtPer filter by fraction
#' @param db data base that was used to profile, default to gtdb.
#' @param candSp candidate species ID
#' @param taxtbl tax table, row.names: 6-digits MIDAS2 ID and `GTDB species` (gtdb)
#' or `Lineage` column.
#' @param filtType "whole" or "group"
#' @param geneType "presabs" or "copynum"
#' @param loadSummary default to FALSE, load summary information.
#' @param loadInfo default to FALSE, load info information.
#' @param loadDepth default to FALSE, load depth information.
#' @param only_stat only samples per species is returned (snpStat and geneStat)
#' @export
#' 
loadMIDAS2 <- function(midas_merge_dir,
                        cl,
                        filtNum=2,
                        db="gtdb",
                        only_stat=FALSE,
                        filtPer=0.8,
                        taxtbl=NULL,
                        candSp=NULL,
                        filtType="group",
                        geneType="copynum",
                        loadSummary=TRUE,
                        loadInfo=TRUE,
                        loadDepth=TRUE) {
  stana <- new("stana")
  if (only_stat) loadSummary <- TRUE
  stana@type <- "MIDAS2"
  stana@db <- db
  stana@cl <- cl

  if (db=="gtdb") {
    tblCol <- "GTDB species"
  } else if (db=="uhgg") {
    tblCol <- "Lineage"
  } else {
    stop("Please provide gtdb or uhgg to `db`")
  }

  if (loadSummary) {
    filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
    snpsSummary <- read.table(filePath, header=1)
    stana@snpsSummary <- snpsSummary
  }

  if (only_stat) {
    filePath <- paste0(midas_merge_dir,"/genes/genes_summary.tsv")
    genesSummary <- read.table(filePath, header=1)

    grnm <- NULL
    for (sn in snpsSummary$sample_name) {
      tmpgr <- NULL
      for (i in names(cl)) {
        if (sn %in% cl[[i]]) {
          tmpgr <- c(tmpgr, i)
        }
      }
      if (length(tmpgr)>1) {
        stop("Sample in multiple groups")
      } else {
        grnm <- c(grnm, tmpgr)
      }
    }

    snpsSummary$group <- grnm

    grnm <- NULL
    for (sn in genesSummary$sample_name) {
      tmpgr <- NULL
      for (i in names(cl)) {
        if (sn %in% cl[[i]]) {
          tmpgr <- c(tmpgr, i)
        }
      }
      if (length(tmpgr)>1) {
        stop("Sample in multiple groups")
      } else {
        grnm <- c(grnm, tmpgr)
      }
    }

    genesSummary$group <- grnm
    snpRet <- snpsSummary |> dplyr::group_by(species_id) |>
        dplyr::count(group)
    geneRet <- genesSummary |> dplyr::group_by(species_id) |>
        dplyr::count(group)

    if (!is.null(taxtbl)) {
      nmdic <- taxtbl[[tblCol]] |> setNames(row.names(taxtbl))
      snpRet$species_name <- nmdic[as.character(snpRet$species_id)]
      geneRet$species_name <- nmdic[as.character(geneRet$species_id)]
    }
    return(
      list("snps"=snpRet,
      "genes"=geneRet)
    )
  }




  stana@mergeDir <- midas_merge_dir
  stana@geneType <- geneType
  stana@sampleFilter <- filtType
  stana@sampleFilterVal <- filtNum
  stana@sampleFilterPer <- filtPer
  snpList <- list()
  snpInfoList <- list()
  snpDepthList <- list()
  geneList <- list()
  clearGn <- NULL
  clearGnSp <- NULL
  clearSn <- NULL
  clearSnSp <- NULL
  snpStat <- NULL
  geneStat <- NULL

  if (filtType=="whole") {
    checkCl <- list("whole"=as.character(unlist(cl)))
  } else if (filtType=="group") {
    checkCl <- cl
  } else {
    stop("Please specify whole or group")
  }
  stana@sampleFilter <- filtType

  qqcat("SNPS\n")
  if (!is.null(candSp)) {specNames <- candSp} else {
    specNames <- list.files(paste0(midas_merge_dir,"/snps"))
  }
  for (i in specNames){
      if (!is.na(as.numeric(i))) {
          qqcat("  @{i}\n")
          if (!is.null(taxtbl)){
              spnm <- taxtbl[i,][[tblCol]]
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
          unlink(paste0(getwd(),"/",cnd))
          ## Info
          if (loadInfo) {
            cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_info.tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d","-f",
                                  paste0(getwd(),"/",cnc),
                                  paste0(getwd(),"/",cnd)),
                    stdout=FALSE, stderr=FALSE)
            info <- read.table(cnd, row.names=1, header=1)
            snpInfoList[[i]] <- info
            unlink(paste0(getwd(),"/",cnd))        
          }
          if (loadDepth) {
            cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_depth.tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d","-f",
                                  paste0(getwd(),"/",cnc),
                                  paste0(getwd(),"/",cnd)),
                    stdout=FALSE, stderr=FALSE)
            depth <- read.table(cnd, row.names=1, header=1)
            snpDepthList[[i]] <- depth
            unlink(paste0(getwd(),"/",cnd))       
          }
          checkPass <- NULL
          spStat <- c(i)
          for (clName in names(checkCl)){
              int <- intersect(colnames(df), checkCl[[clName]])
              spStat <- c(spStat, int |> length())
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(checkCl[[clName]])*filtPer)
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
          snpStat <- rbind(snpStat, spStat)
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
              spnm <- taxtbl[i,][[tblCol]]
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
          spStat <- c(i)
          for (clName in names(checkCl)){
              int <- intersect(colnames(df), checkCl[[clName]])
              spStat <- c(spStat, int |> length())
              qqcat("      Number of samples in @{clName}: @{length(int)}\n")
              checkPass <- c(checkPass, length(int)>filtNum | length(int)>length(checkCl[[clName]])*filtPer)
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
          geneStat <- rbind(geneStat, spStat)
          unlink(paste0(getwd(),"/",cnd))
      }
  }
  snpStat <- snpStat |> data.frame() |> `colnames<-`(c("species",names(checkCl)))
  row.names(snpStat) <- snpStat$species
  for (i in names(checkCl)) {
    snpStat[[i]] <- as.numeric(snpStat[[i]])
  }
  geneStat <- geneStat |> data.frame() |> `colnames<-`(c("species",names(checkCl)))
  row.names(geneStat) <- geneStat$species
  for (i in names(checkCl)) {
    geneStat[[i]] <- as.numeric(geneStat[[i]])
  }

  stana@freqTableSnps <- snpStat
  stana@freqTableGenes <- geneStat
  if (!is.null(clearSn)) stana@clearSnps <- clearSn
  if (!is.null(clearGn)) stana@clearGenes <- clearGn
  stana@snps <- snpList
  stana@snpsInfo <- snpInfoList
  stana@snpsDepth <- snpDepthList
  stana@genes <- geneList
  stana@ids <- union(names(geneList),names(snpList))
  stana <- initializeStana(stana,cl)
  if (!is.null(taxtbl)){
    if (!is.null(clearSnSp)) stana@clearSnpsSpecies <- clearSnSp
    if (!is.null(clearGnSp)) stana@clearGenesSpecies <- clearGnSp
  }
  stana
}