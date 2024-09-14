#' loadMIDAS2
#' 
#' load the MIDAS2 merge command output.
#' The location to lz4 binary must be added to PATH.
#' Taxonomy table can be loaded from downloaded MIDAS2 db directory (metadata.tsv).
#' If provided, additionally show tax names.
#' 
#' @param midas_merge_dir path to merged directory
#' @param cl named list of category for samples
#' @param filtBy filter by {snps} number or {genes} number (default to snps)
#' @param filtNum the species with the samples above this number will be returned
#' @param filtPer filter by fraction
#' @param db data base that was used to profile, default to gtdb.
#' @param candSp candidate species ID
#' @param filtType "whole" or "group"
#' @param geneType "presabs" or "copynum"
#' @param loadSummary default to FALSE, load summary information.
#' @param loadInfo default to FALSE, load info information.
#' @param loadDepth default to FALSE, load depth information.
#' @param only_stat only samples per species is returned (snpStat and geneStat)
#' @importFrom dplyr mutate
#' @export
#' 
loadMIDAS2 <- function(midas_merge_dir,
                        cl=NULL,
                        filtNum=2,
                        db="gtdb",
                        only_stat=FALSE,
                        filtBy="snps",
                        filtPer=0.8,
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
    if (!is.null(cl)) {
        stana@cl <- cl
    }

    if (db=="gtdb") {
        tblCol <- "GTDB species"
    } else if (db=="uhgg") {
        tblCol <- "Lineage"
    } else {
        stop("Please provide gtdb or uhgg to `db`")
    }

    ## Load summary for filtering
    ## Probably filter based on these summaries, but what if 
    ## These files are not available and only species directories are present?
    filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
    snpsSummary <- read.table(filePath, header=1)
    stana@snpsSummary <- snpsSummary
    filePath <- paste0(midas_merge_dir,"/genes/genes_summary.tsv")
    genesSummary <- read.table(filePath, header=1)
    stana@genesSummary <- genesSummary

    if (is.null(cl)) {
        cl <- list("no_group"=unique(c(genesSummary$sample_name, snpsSummary$sample_name)))
        stana@cl <- cl
    }

    changer <- listToNV(cl)
    ## error in duplicate
    snpsSummary$group <- changer[snpsSummary$sample_name]
    genesSummary$group <- changer[genesSummary$sample_name]
    snpRet <- snpsSummary |> dplyr::group_by(species_id) |>
        dplyr::count(group)
    geneRet <- genesSummary |> dplyr::group_by(species_id) |>
        dplyr::count(group)
    ## Change tax name if available
    
    snpRet[["species_name"]] <- loadDic()[[db]][as.character(snpRet$species_id)] 
    geneRet[["species_name"]] <- loadDic()[[db]][as.character(geneRet$species_id)]
    
    snpRet$species_id <- as.character(snpRet$species_id)
    geneRet$species_id <- as.character(geneRet$species_id)

    if (only_stat) {
    	dics <- loadDic()
        return(
            list("snps"=snpRet %>% mutate(species_description=dics[[db]][species_id]),
                "genes"=geneRet %>% mutate(species_description=dics[[db]][species_id]))
        )        
    }

    stana@mergeDir <- midas_merge_dir
    stana@geneType <- geneType
    stana@sampleFilter <- filtType
    stana@sampleFilterVal <- filtNum
    stana@sampleFilterPer <- filtPer

    clearGn <- NULL
    clearGnSp <- NULL
    clearSn <- NULL
    clearSnSp <- NULL

    snpStat <- NULL
    geneStat <- NULL

    if (!is.null(candSp)) {
        ispecNames <- candSp
    } else {
        if (filtBy=="snps") {
            checker <- snpRet
        } else {
            checker <- geneRet
        }
        if (filtType=="whole") {
            wn <- length(unlist(cl))
            ispecNames <- checker %>% group_by(species_id) %>%
                summarise(grn=sum(n)) %>%
                mutate(fnum=grn>filtNum, fper=grn>wn*filtPer) %>%
                filter(fnum & fper) %>% pull(species_id) %>% unique()
        } else if (filtType=="group") {
            ispecNames <- lapply(names(cl), function(cn) {
                gn <- length(cl[[cn]])
                checker %>% 
                    filter(.data$group==cn) %>%
                    mutate(fnum=n>filtNum, fper=n>gn*filtPer) %>%
                    filter(fnum & fper) %>% pull(species_id) %>% unique()
            })
            ispecNames <- Reduce(intersect, ispecNames)
        } else {
            stop("Please specify whole or group")
        }
        stana@sampleFilter <- filtType
    }

    ## Ensure that the directory is present
    snpsSpecNames <- list.files(paste0(midas_merge_dir,"/snps"))
    specNames <- intersect(snpsSpecNames, ispecNames)
    if (length(specNames)!=0) stana@clearSnps <- specNames

    snpsList <- lapply(specNames, function(i) {
        if (!is.na(as.numeric(i))) {## Check that ID is number
            qqcat("  @{i}\n")
            spnm <- loadDic()[[db]][i]
            qqcat("  @{spnm}\n")
            cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_freqs.tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d","-f",
                                paste0(getwd(),"/",cnc),
                                paste0(getwd(),"/",cnd)),
                    stdout=FALSE, stderr=FALSE)
            # df <- read.table(cnd, row.names=1, header=1)
            df <- data.frame(fread(cnd))
            row.names(df) <- df$site_id
            df$site_id <- NULL
            qqcat("    Number of snps: @{dim(df)[1]}\n")
            qqcat("    Number of samples: @{dim(df)[2]}\n")
            unlink(paste0(getwd(),"/",cnd))
            return(df)
        }
    }) %>% setNames(specNames)

    ## Info
    if (loadInfo) {
        snpInfoList <- lapply(specNames, function(i) {
            cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_info.tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d","-f",
                            paste0(getwd(),"/",cnc),
                            paste0(getwd(),"/",cnd)),
                stdout=FALSE, stderr=FALSE)
            info <- read.table(cnd, row.names=1, header=1)
            unlink(paste0(getwd(),"/",cnd))
            return(info)
        }) %>% setNames(specNames)
        stana@snpsInfo <- snpInfoList
    }
     
    if (loadDepth) {
        snpDepthList <- lapply(specNames, function(i) {
            cnc <- paste0(midas_merge_dir,"/snps/",i,"/",i,".snps_depth.tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d","-f",
                    paste0(getwd(),"/",cnc),
                    paste0(getwd(),"/",cnd)),
                    stdout=FALSE, stderr=FALSE)
            depth <- read.table(cnd, row.names=1, header=1)
            unlink(paste0(getwd(),"/",cnd))
            return(depth)            
        }) %>% setNames(specNames)
        stana@snpsDepth <- snpDepthList  
    }

    ## Ensure that the directory is present
    genesSpecNames <- list.files(paste0(midas_merge_dir,"/genes"))
    specNames <- intersect(genesSpecNames, ispecNames)
    if (length(specNames)!=0) stana@clearGenes <- specNames

    geneList <- lapply(specNames, function(i) {
        if (!is.na(as.numeric(i))) {
            qqcat("  @{i}\n")
            spnm <- loadDic()[[db]][i]
            qqcat("  @{spnm}\n")
            cnc <- paste0(midas_merge_dir,
                        "/genes/",i,"/",i,".genes_",geneType,".tsv.lz4")
            cnd <- gsub(".lz4","",cnc)
            system2("lz4", args=c("-d", "-f",
                                paste0(getwd(),"/",cnc),
                                paste0(getwd(),"/",cnd)),
                    stdout=FALSE, stderr=FALSE)
            df <- data.frame(fread(cnd))
            row.names(df) <- df$gene_id
            df$gene_id <- NULL
            qqcat("    Number of genes: @{dim(df)[1]}\n")
            qqcat("    Number of samples: @{dim(df)[2]}\n")
            unlink(paste0(getwd(),"/",cnd))
            df
        }
    }) %>% setNames(specNames)



    stana@snps <- snpsList
    stana@genes <- geneList
    stana@ids <- union(names(geneList),names(snpsList))
    stana@names <- loadDic()[[db]][stana@ids] %>% setNames(stana@ids)
    stana <- initializeStana(stana,cl)


    stana@clearSnpsSpecies <- loadDic()[[db]][stana@clearSnps]
    stana@clearGenesSpecies <- loadDic()[[db]][stana@clearGenes]

    stana
}