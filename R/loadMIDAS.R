
#' loadMIDAS
#' 
#' Assess and store profile for species and return filtered species 
#' based on the number of samples for each category or whole population.
#' For MIDAS only.
#'
#' @param midas_merge_dir output directory of merge_midas.py
#' @param cl named list of sample IDs
#' @param filtType "whole" or "group"
#' @param geneType "presabs" or "copynum"
#' @param filtNum The species with number above this threshold
#'                for each category is returned
#' @param filtPer filter by fraction
#' @param candSp candidate species ID
#' @param loadSummary default to FALSE, load summary information.
#' @param loadInfo default to FALSE, load info information.
#' @param loadDepth default to FALSE, load depth information.
#' @param only_stat return stat for snp and gene only
#' @import GetoptLong
#' @export
loadMIDAS <- function(midas_merge_dir,
  cl, filtType="group", candSp=NULL,
  only_stat=FALSE,
  filtNum=2, filtPer=0.8,
  loadSummary=TRUE,
  loadDepth=TRUE,
  loadInfo=TRUE,
  geneType="copynum") {
  stana <- new("stana")
  stana@type <- "MIDAS1"
  stana@mergeDir <- midas_merge_dir
  stana@sampleFilter <- filtType
  stana@sampleFilterVal <- filtNum
  stana@sampleFilterPer <- filtPer
  stana@cl <- cl
  dirLs <- list.files(midas_merge_dir)
  specNames <- NULL
  for (d in dirLs) {
    if (dir.exists(paste0(midas_merge_dir,"/",d))){
      specNames <- c(specNames, d)
    }
  }
  if (only_stat) {
  		statssnps <- NULL
		statsgenes <- NULL
  		for (sp in specNames) {
			ret <- NULL;
			ret2 <- NULL;
  			if ("snps_summary.txt" %in% list.files(paste0(midas_merge_dir,"/",sp))){
				snp_df <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_summary.txt"),
					row.names=1, sep="\t", header=1)
					ret <- c(ret, sp)
					snpsmp <- row.names(snp_df)
					for (clnm in names(cl)) {
						ret <- c(ret, sum(cl[[clnm]] %in% snpsmp))
					}
  			}
  			if ("genes_summary.txt" %in% list.files(paste0(midas_merge_dir,"/",sp))){
				gene_df <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_summary.txt"),
				row.names=1, sep="\t", header=1)
				genessmp <- row.names(gene_df)
				ret2 <- c(ret2, sp)
				for (clnm in names(cl)) {
					ret2 <- c(ret2, sum(cl[[clnm]] %in% genessmp))
				}
			}
			statssnps <- rbind(statssnps, ret)
			statsgenes <- rbind(statsgenes, ret2)
  		}
  		snpsdf <- rbind(statssnps) |> data.frame() |> `colnames<-`(c("species", names(cl)))
  		row.names(snpsdf) <- snpsdf$species
  		genesdf <- rbind(statsgenes) |> data.frame() |> `colnames<-`(c("species", names(cl)))
  		row.names(genesdf) <- genesdf$species
  		return(list(
  			"snps"=snpsdf,
  			"genes"=genesdf
  			))

  }

  stana@geneType <- geneType
  clearSn <- NULL
  clearGn <- NULL
  freqtblSn <- NULL
  freqtblGn <- NULL
  snpList <- list()
  geneList <- list()
  ## Summary is loaded per species and concatenated
  snpInfoList <- list()
  snpDepthList <- list()
  snpSummary <- list()
  if (filtType=="whole") {
    checkCl <- list("whole"=as.character(unlist(cl)))
  } else if (filtType=="group") {
    checkCl <- cl
  } else {
    stop("Please specify whole or group")
  }
  stana@sampleFilter=filtType

  if (!is.null(candSp)) {specNames <- candSp}
  stana@ids <- specNames
  for (sp in specNames){
    pnum <- c(sp)
    qqcat("@{sp}\n")
    qqcat("  Snps\n")
    cont <- list.files(paste0(midas_merge_dir,"/",sp))
    if (loadSummary & "snps_summary.txt" %in% cont){
      summ <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_summary.txt"),
                         sep="\t",header=1)
      summ$species_id <- sp
      snpSummary[[sp]] <- summ

    }
    if (loadInfo & "snps_info.txt" %in% cont) {
      info <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_info.txt"),
                         sep="\t",header=1,row.names=1)
      ## Rename the ID row.names
      row.names(info) <- paste0(info$ref_id, "_", info$ref_pos)
      snpInfoList[[sp]] <- info

    }
    if (loadDepth & "snps_depth.txt" %in% cont) {
      depth <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_depth.txt"),
                         sep="\t",header=1,row.names=1)
      snpDepthList[[sp]] <- depth   
    }


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
  if (!is.null(clearSn)) stana@clearSnps <- clearSn
  if (!is.null(clearGn)) stana@clearGenes <- clearGn
  stana@freqTableSnps <- freqtblSn
  stana@freqTableGenes <- freqtblGn
  stana@snps <- snpList
  stana@genes <- geneList
  stana@snpsInfo <- snpInfoList
  stana@snpsDepth <- snpDepthList
  stana@snpsSummary <- do.call(rbind, snpSummary)

  stana <- initializeStana(stana,cl)

  qqcat("Overall, @{length(clearSn)} species met criteria in SNPs\n")
  qqcat("Overall, @{length(clearGn)} species met criteria in genes\n")
  stana
}
