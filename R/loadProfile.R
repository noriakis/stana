
#' loadInStrain
#' 
#' For InStrain, you should provide `compare` output
#' with --bams options implmented from version 1.6.
#' It computes MAF of var_base per position and return the frequency matrix.
#' Note it considers con_base and var_base in info table for calculation.
#' 
#' @param compare_out_dir output directory of compare.
#' need `output` in the directory
#' @param candidate_species candidate species ID, e.g. `GUT_GENOME000024`
#' @param cl named list of grouping
#' @param just_species if TRUE, returns the vector of species IDs
#' @param fill_na fill NA in resulting data.frame by -1
#' 
#' @export
loadInStrain <- function(compare_out_dir,
                         candidate_species,
                         cl=NULL,
                         just_species=FALSE,
                         fill_na=TRUE) {
    sample_threshold <- 1
    output_list <- list.files(paste0(compare_out_dir,"/output"))
    stana <- new("stana")
    if (!is.null(cl)) {stana@cl <- cl}
    stana@type <- "InStrain"
    stana@mergeDir <- compare_out_dir
    snps <- list()
    snpsInfoList <- list()

    ctpath <- paste0(compare_out_dir,"/output/",output_list[grepl("comparisonsTable.tsv", output_list)])
    ct <- data.table::fread(ctpath)
    stana@comparisonTable <- ct
    gcpath <- paste0(compare_out_dir,"/output/",output_list[grepl("genomeWide_compare.tsv", output_list)])
    gc <- data.table::fread(gcpath)
    stana@genomeWideCompare <- gc
    scpath <- paste0(compare_out_dir,"/output/",output_list[grepl("strain_clusters.tsv", output_list)])
    sc <- data.table::fread(scpath)
    stana@strainClusters <- sc

    if (sum(grepl("pooled_SNV_data",output_list))!=0) {
        if (!just_species) {
          qqcat("Loading allele count table\n")
          path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data.tsv.gz",output_list)])
          tbl <- data.table::fread(path)          
        }

        keys_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data_keys.tsv",output_list)])
        qqcat("Loading key table\n")
        keys <- data.table::fread(keys_path)
        info_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_info.tsv.gz",output_list)])
        qqcat("Loading info table\n")
        info <- data.table::fread(info_path)

        info$scaffold_position <- paste0(info$scaffold,"|",info$position)


        varBase <- info$var_base
        names(varBase) <- info$scaffold_position
        conBase <- info$con_base
        names(conBase) <- info$scaffold_position
        
        sps <- unique(paste0(sapply(strsplit(keys$scaffold, "_"),"[",1),"_",sapply(strsplit(keys$scaffold, "_"),"[",2)))
        if (just_species) {
            return(sps)
        }

        snpsInfoList[[candidate_species]] <- info[grepl(candidate_species, info$scaffold),]
        stana@ids <- sps
        qqcat("Candidate species: @{candidate_species}\n")
        candKeys <- keys[grepl(candidate_species,keys$scaffold),]$key
        qqcat("  Candidate key numbers: @{length(candKeys)}\n")
        ret <- function(x) {
            smp <- x[1]
            scaff <- x[2]
            scp <- paste0(scaff,"|",as.numeric(x[3]))
            det <- c(smp, scp)
            cb <- x[8]
            vb <- x[9]
            
            read_num <- as.numeric(x[4:7])
            names(read_num) <- names(x[4:7])
            ## Return -1 in zero depth
            if (sum(read_num)==0) {val <- -1} else {
                val <- as.numeric((read_num / sum(read_num))[vb])
            }
            det <- c(det, cb, vb, val)
            det
        }
        mafs <- NULL
        
        sdt <- tbl[tbl$scaffold %in% candKeys,]
        qqcat("  Dimension of pooled SNV table for species: @{dim(sdt)[1]}\n")
        # Only bi-allelic position
        # sdt <- sdt[apply(sdt[,4:7],1,function(x)sum(x==0))==2,]

        nm_sample <- keys$sample
        names(nm_sample) <- keys$key
        nm_sample <- nm_sample[nm_sample!=""]
        nm_scaffold <- keys$scaffold
        names(nm_scaffold) <- keys$key
        nm_scaffold <- nm_scaffold[nm_scaffold!=""]
        
        sdt$sample <- nm_sample[as.character(sdt$sample)]
        sdt$scaffold <- nm_scaffold[as.character(sdt$scaffold)]
        sdt$conBase <- conBase[paste0(sdt$scaffold,"|",sdt$position)]
        sdt$varBase <- varBase[paste0(sdt$scaffold,"|",sdt$position)]
        

        qqcat("Calculating MAF\n")
        maf <- apply(sdt, 1, function(x) ret(x))
        maf <- data.frame(t(maf))
        
        maf <- maf |> `colnames<-`(c("id","scaffold_position","major_allele","minor_allele","maf"))
        maf$maf <- as.numeric(maf$maf)
        maf$scaffold_position <- paste0(maf$scaffold_position,"|",maf$major_allele,">",maf$minor_allele)

        sample_snv <- tidyr::pivot_wider(maf, id_cols=id,
                                 names_from = scaffold_position,
                                 values_from = maf)
        thresh <- apply(sample_snv[2:ncol(sample_snv)], 2, function(x) sum(!is.na(x))) > sample_threshold
        snvDf <- data.frame(sample_snv[,c("id",names(thresh[thresh]))], check.names=FALSE)
        row.names(snvDf) <- snvDf$id
        snvDf$id <- NULL
        snvDf <- data.frame(t(snvDf), check.names=FALSE)
        ## Return the data.frame
        if (fill_na) {
          snvDf[is.na(snvDf)] <- -1
        }
        snps[[candidate_species]] <- snvDf
        stana@snps <- snps
        stana@snpsInfo <- snpsInfoList
        stana <- initializeStana(stana,cl)
    } else {
        qqcat("No pooled results present\n")
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
#' @param filtNum the species with the samples above this number will be returned
#' @param filtPer filter by fraction
#' @param filtType "whole" or "group"
#' @param candSp candidate species ID
#' @import GetoptLong
#' @export
loadmetaSNV <- function(metasnv_out_dir, cl=NULL,
                        filtType="group", candSp=NULL,
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
#' @param geneType "presabs" or "copynum"
#' @param filtNum The species with number above this threshold
#'                for each category is returned
#' @param filtPer filter by fraction
#' @param candSp candidate species ID
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
#' @export
#' 
loadMIDAS2 <- function(midas_merge_dir,
                        cl,
                        filtNum=2,
                        db="gtdb",
                        filtPer=0.8,
                        taxtbl=NULL,
                        candSp=NULL,
                        filtType="group",
                        geneType="copynum",
                        loadSummary=FALSE,
                        loadInfo=FALSE,
                        loadDepth=FALSE) {
  stana <- new("stana")
  stana@type <- "MIDAS2"
  stana@db <- db
  stana@cl <- cl
  if (loadSummary) {
    filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
    snpsSummary <- read.table(filePath, header=1)
    stana@snpsSummary <- snpsSummary
  }
  if (db=="gtdb") {
    tblCol <- "GTDB species"
  } else if (db=="uhgg") {
    tblCol <- "Lineage"
  } else {
    stop("Please provide gtdb or uhgg to `db`")
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
  stana@clearSnps <- clearSn
  stana@clearGenes <- clearGn
  stana@snps <- snpList
  stana@snpsInfo <- snpInfoList
  stana@snpsDepth <- snpDepthList
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
