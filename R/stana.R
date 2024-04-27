#' Note that gene summary and depth is not loaded in MIDAS and MIDAS2.
#' @import methods
#' @importFrom utils object.size
#' @importFrom grDevices colorRampPalette
#' @slot ids identifiers distinguishing species
#' @slot names name of species (named vector)
#' @slot snpsInfo snp info
#' @slot snpsDepth snp depth
#' @slot snpsSummary snp summary
#' @slot kos stores ko table
#' @slot type type of pipeline, such as `MIDAS`
#' @slot eggNOG list of path for eggNOG mapper v2 results for species
#' @slot map slot storing the mapping data.frame for gene ID and orthology
#' @slot coefMat slot storing the strain (or subspecies) abundances
#' @slot includeSNVID if SNV allele frequency is used in the calculation,
#' the IDs in the slot is subset by default.
setClass("stana", slots=list(
                            type="character",
                            cl="list",
                            meta="data.frame",
                            mergeDir="character",
                            db="character",
                            ids="character",
                            names="character",
                            comparisonTable="data.table",
                            genomeWideCompare="data.table",
                            strainClusters="data.table",
                            snps="list",
                            eggNOG="list",
                            NMF = "list",
                            kos="list",
                            includeSNVID="list",
                            snpsInfo="list",
                            snpsDepth="list",
                            snpsSummary="data.frame",
                            genesSummary="data.frame",
                            relab="data.frame",
                            geneType="character",
                            genes="list",
                            clearSnps="vector",
                            clearGenes="vector",
                            clearSnpsSpecies="vector",
                            clearGenesSpecies="vector",
                            freqTableSnps="data.frame",
                            freqTableGenes="data.frame",
                            fastaList="list",
                            treeList="list",
                            adonisList="list",
                            treePlotList="list",
                            sampleFilter="character",
                            sampleFilterVal="numeric",
                            sampleFilterPer="numeric",
                            heatmap="list",
                            coefMat="list",
                            colors="vector",
                            geneCluster="list",
                            paFilterUp="numeric",
                            paFilterDown="numeric",
                            gsea="list",
                            map="list"))

#' cat_subtle
#' @noRd
cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))

#' show
#' print the description of stana
#' @importFrom dplyr group_by mutate summarise n
#' @noRd
setMethod("show",
  signature(object="stana"),
  function(object) {
  	cat_subtle("# A stana: ", object@type, "\n", sep="")
  	if (!identical(object@db, character(0))) {
  		cat_subtle("# Database: ", object@db, "\n", sep="")
  	}
    cat_subtle("# Loaded directory: ", object@mergeDir, "\n", sep="")
    cat_subtle("# Species number: ", length(object@ids), "\n", sep="")
    if (length(object@cl)!=0) {
        cat_subtle("# Group info (list): ", paste0(names(object@cl), collapse="/"), "\n", sep="")	
    }
    if (dim(object@meta)[1]!=0) {
        cat_subtle("# Group column (DF): ", paste0(colnames(object@meta), collapse="/"), "\n", sep="")	
    }
    # if (object@type %in% c("MIDAS","MIDAS2")) {
    #   qqcat("Filter type: @{object@sampleFilter}, number: @{object@sampleFilterVal}, proportion: @{object@sampleFilterPer}\n")
    # }
    if (length(object@snps)!=0) {
      cat_subtle("# Loaded SNV table: ", length(object@snps), " ID: ", paste0(names(object@snps)[1], collapse="/"), "\n", sep="")
    }
    # if (object@type %in% c("MIDAS","MIDAS2")) {
    #   qqcat("  Species cleared SNV filter: @{length(object@clearSnps)}\n")
    # }
    if (length(object@genes)!=0) {
      cat_subtle("# Loaded gene table: ", length(object@genes), " ID: ", paste0(names(object@genes)[1], collapse="/"), "\n", sep="")
    }
    # if (object@type %in% c("MIDAS","MIDAS2")) {
    #   qqcat("  Species cleared gene filter: @{length(object@clearGenes)}\n")
    # }
    if (length(object@kos)!=0) {
      cat_subtle("# Loaded KO table: ", length(object@kos), " ID: ", paste0(names(object@kos)[1], collapse="/"), "\n", sep="")
    }
    if (length(object@fastaList)!=0) {
      cat_subtle("# Inferred fasta: ", length(object@kos), " ID: ", paste0(names(object@kos)[1], collapse="/"), "\n", sep="")
    }
    cat_subtle("# Size: ", object.size(object), " B\n", sep="")
  })

#' loadDic
#' @noRd
loadDic <- function() {
    uhgg <- readRDS(system.file("extdata", "uhggdic.rds", package = "stana"))
    gtdb <- readRDS(system.file("extdata", "gtdbdic.rds", package = "stana"))
    list("uhgg"=uhgg, "gtdb"=gtdb)
}


#' summary
#' print summary information
#' @param object stana object
#' @param ... other arguments
#' @export
setGeneric("summary", function(object, ...) standardGeneric("summary"))
#' summary
#' print summary information
#' @param object stana object
#' @param ... other arguments
#' @export
setMethod("summary", "stana",
    function(object, ...) {
        if (object@type %in% c("MIDAS", "MIDAS2")) {
            cat_subtle("# \n")
            cat_subtle("# SNV description\n")
            df <- object@snpsSummary %>%
                dplyr::filter(.data$species_id %in% names(object@snps))
            if (object@type=="MIDAS2") {
                df$species_id <- loadDic()[[object@db]][as.character(df$species_id)]
            }
            if (length(object@cl)!=0) {
                print(df %>% mutate(group=listToNV(object@cl)[sample_name]) %>%
                group_by(group, species_id) %>%
                summarise(n=n()))
            } else {
                print(df %>%
                group_by(species_id) %>%
                summarise(n=n()))
            }
            cat_subtle("# Gene description\n")
            df <- object@genesSummary %>%
                dplyr::filter(.data$species_id %in% names(object@genes))
            if (object@type=="MIDAS2") {
                df$species_id <- loadDic()[[object@db]][as.character(df$species_id)]
            }
            if (length(object@cl)!=0) {
                print(df %>% mutate(group=listToNV(object@cl)[sample_name]) %>%
                group_by(group, species_id) %>%
                summarise(n=n()))
            } else {
                print(df %>%
                group_by(species_id) %>%
                summarise(n=n()))
            }
        }
    }    
)



#' siteFilter
#' @param x stana object
#' @param sp species ID
#' @param exp expression for filtering the snps info
#' @export
setGeneric("siteFilter", function(x, sp, exp) standardGeneric("siteFilter"))

#' siteFilter
#' @param x stana object
#' @param sp species ID
#' @param exp expression for filtering the snps info
#' @export
setMethod("siteFilter", "stana",
    function(x, sp, exp) {
    	info <- x@snpsInfo[[sp]]
    	ret <- info %>% 
    	    dplyr::filter(!!enquo(exp)) %>%
    	    row.names(.)
    	cat_subtle("# total of ", length(ret), " obtained from ", dim(info)[1], "\n", sep="")
    	x@includeSNVID[[sp]] <- ret
    	x
})


#' setSNVID
#' @param x stana object
#' @param sp species ID
#' @param IDs SNV ID
#' @export
setGeneric("setSNVID", function(x, sp, IDs) standardGeneric("setSNVID"))

#' setSNVID
#' @param x stana object
#' @param sp species ID
#' @param IDs SNV ID
#' @export
setMethod("setSNVID", "stana",
    function(x, sp, IDs) {
    x@includeSNVID[[sp]] <- IDs
    x
})



#' check
#' check and output statistics based on conditional formulas for summary
#' @param x stana boject
#' @param exp expression for filtering the snps summary
#' @param target snps or genes
#' @export
setGeneric("check", function(x, exp, target="snps") standardGeneric("check"))
#' check
#' check and output statistics based on conditional formulas for summary
#' @param x stana boject
#' @param exp expression for filtering the snps summary
#' @param target snps or genes
#' @export
setMethod("check", "stana",
    function(x, exp, target) {
    	if (target=="snps") {
        	if (dim(x@snpsSummary)[1]==0) {
        		stop("Please provide SNPs summary for filtering")
        	}
        	ret <- x@snpsSummary %>% 
        	    dplyr::filter(!!enquo(exp)) %>% dplyr::group_by(species_id) %>%
        	    dplyr::summarize(n=dplyr::n())
    	}
    	if (target=="genes") {
        	if (dim(x@genesSummary)[1]==0) {
        		stop("Please provide genes summary for filtering")
        	}
        	ret <- x@genesSummary %>% 
        	    dplyr::filter(!!enquo(exp)) %>% dplyr::group_by(species_id) %>%
        	    dplyr::summarize(n=dplyr::n())
    	}
    	if (x@type=="MIDAS2") {
    		ret <- ret %>% dplyr::mutate(species_description=loadDic()[[x@db]][as.character(species_id)])    		
    	}
    	return(ret)
    })


#' filter
#' filter the stana object based on species ID
#' @param x stana object
#' @param ids species ID
#' @param target snps or genes
#' @export
setGeneric("filter", function(x, ids, target="snps") standardGeneric("filter"))

#' filter
#' filter the stana object based on species ID
#' @param x stana object
#' @param ids species ID
#' @param target snps or genes
#' @export
setMethod("filter", "stana",
    function(x, ids, target) {
        ids <- intersect(x@ids, ids)
    	if (target=="snps") {
            ls <- x@snps[ids]
            ls <- ls[lapply(ls, function(x) !is.null(x)) %>% unlist()]
    		x@snps <- ls
            x@ids <- ids
            x@names <- x@names[ids]
    	}
    	if (target=="genes") {
            ls <- x@genes[ids]
            ls <- ls[lapply(ls, function(x) !is.null(x)) %>% unlist()]
            x@genes <- ls
            x@ids <- ids
            x@names <- x@names[ids]
    	}
    	return(x)
    })


#' getSlot 
#' get the slot values
#' @param x stana object
#' @param slot slot name
#' @export
setGeneric("getSlot",
    function(x, slot) standardGeneric("getSlot"))

#' getSlot 
#' get the slot values
#' @param x stana object
#' @param slot slot name
#' @export
setMethod("getSlot", "stana",
    function(x, slot) attr(x, slot))

#' getFasta
#' get fasta list from stana
#' @param x stana object
#' @export
setGeneric("getFasta",
    function(x) standardGeneric("getFasta"))

#' getFasta
#' get fasta list from stana
#' @param x stana object
#' @export
setMethod("getFasta", "stana",
    function(x) attr(x, "fastaList"))

#' getTree
#' get tree list from stana
#' @param x stana object
#' @export
setGeneric("getTree",
    function(x) standardGeneric("getTree"))

#' getTree
#' get tree list from stana
#' @param x stana object
#' @export
setMethod("getTree", "stana",
    function(x) attr(x, "treeList"))

#' getTreePlot
#' get tree plot list from stana
#' @param x stana object
#' @export
setGeneric("getTreePlot",
    function(x) standardGeneric("getTreePlot"))

#' getTreePlot
#' get tree plot list from stana
#' @param x stana object
#' @export
setMethod("getTreePlot", "stana",
    function(x) attr(x, "treePlotList"))

#' getID
#' get ID list from stana
#' @param x stana object
#' @export
setGeneric("getID",
    function(x) standardGeneric("getID"))


#' getID
#' get ID list from stana
#' @param x stana object
#' @export
setMethod("getID", "stana",
    function(x) attr(x, "ids"))

#' getAdonis
#' get adonis list from stana
#' @param x stana object
#' @export
setGeneric("getAdonis",
    function(x) standardGeneric("getAdonis"))

#' getAdonis
#' get adonis list from stana
#' @param x stana object
#' @export
setMethod("getAdonis", "stana",
    function(x) attr(x, "adonisList"))

#' getCl
#' get grouping information list from stana
#' @param x stana object
#' @export
setGeneric("getCl",
    function(x) standardGeneric("getCl"))

#' getCl
#' get grouping information list from stana
#' @param x stana object
#' @export
setMethod("getCl", "stana",
    function(x) attr(x, "cl"))

#' getGeneID
#' get gene IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setGeneric("getGeneID",
    function(x, candSp) standardGeneric("getGeneID"))

#' getGeneID
#' get gene IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setMethod("getGeneID", "stana",
    function(x, candSp) {
        if (is.null(stana@genes[[candSp]])) {
            stop("No gene table available")
        } else {
            return(row.names(stana@genes[[candSp]]))
        }
})

#' getSNVID
#' get SNV IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setGeneric("getSNVID",
    function(x, candSp) standardGeneric("getSNVID"))

#' getSNVID
#' get SNV IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setMethod("getSNVID", "stana",
    function(x, candSp) {
        if (is.null(stana@snps[[candSp]])) {
            stop("No SNV table available")
        } else {
            return(row.names(stana@snps[[candSp]]))
        }
})

#' getKOID
#' get gene family IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setGeneric("getKOID",
    function(x, candSp) standardGeneric("getKOID"))

#' getKOID
#' get gene family IDs from stana
#' @param x stana object
#' @param candSp species
#' @export
setMethod("getKOID", "stana",
    function(x, candSp) {
        if (is.null(stana@kos[[candSp]])) {
            stop("No KO table available")
        } else {
            return(row.names(stana@kos[[candSp]]))
        }
})



#' setSlot
#' general method for setting the slot
#' @param x stana object
#' @param y slot name
#' @param z set variable
#' @export
setGeneric("setSlot",
    function(x, y, z) standardGeneric("setSlot"))

#' setSlot
#' general method for setting the slot
#' @param x stana object
#' @param y slot name
#' @param z set variable
#' @export
setMethod("setSlot", "stana",
    function(x, y, z) {
        attr(x, y) <- z
        x
})



#' initializeStana
#' @noRd
initializeStana <- function(stana,cl) {
  stana@colors <- getColors(cl)
  # faList <- vector("list", length(stana@ids))
  # names(faList) <- stana@ids
  # stana@fastaList <- faList
  # treeList <- vector("list", length(stana@ids))
  # names(treeList) <- stana@ids
  # stana@treeList <- treeList
  # treePlotList <- vector("list", length(stana@ids))
  # names(treePlotList) <- stana@ids
  # stana@treePlotList <- treePlotList
  # adonisList <- vector("list", length(stana@ids))
  # names(adonisList) <- stana@ids
  # stana@adonisList <- adonisList
  stana@treeList <- list()
  stana@treePlotList <- list()
  stana@fastaList <- list()
  stana@adonisList <- list()
  stana
}

#' setAnnotation
#' @param stana stana object
#' @param annotList named list of the path to eggNOG annotation file
#' (*emapper.annotations)
#' @export
#' @return stana object
setAnnotation <- function(stana, annotList) {
	stana@eggNOG <- annotList
	return(stana)
}


#' setMap
#' @param stana stana object
#' @param candSp species ID
#' @param map mapping data.frame, two column consisting of ID and value.
#' The function internally changes the column name.
#' ID is gene ID in `genes` slot, and value indicates orthology ID (like K00001)
#' @export
#' @return stana object
setMap <- function(stana, candSp, map) {
	if (dim(map)[2]!=2) {stop("Please provide two-column data.frame")}
	colnames(map) <- c("ID", "value")
	stana@map[[candSp]] <- map
	return(stana)
}


#' getColors
#' @import RColorBrewer
#' @noRd
getColors <- function(cl){
  numgr <- length(names(cl))
  if (numgr > 2) {
    cols <- brewer.pal(numgr, "PuOr") 
  } else if (numgr == 2) {
    three <- brewer.pal(3, "PuOr")
    cols <- c(three[1], three[3])
  } else {
    cols <- brewer.pal(3, "PuOr")[1]
  }
}

#' changeColors
#' @param stana stana object
#' @param colors color vector
#' @export
#' @return stana object
changeColors <- function(stana, colors) {
	stana@colors <- colors
	return(stana)
}

#' snv
#' @param dfs named list of dfs
#' @export
snv <- function(dfs) {
    stana <- new("stana")
    stana@type <- "manual"

    if (is.null(names(dfs))) {
        ids <- paste0("species", seq_len(length(dfs)))
    } else {
        ids <- names(dfs)
    }
    stana@snps <- dfs
    stana@ids <- ids
    return(stana)
}

#' gene
#' @param dfs named list of dfs
#' @export
gene <- function(dfs) {
    stana <- new("stana")
    stana@type <- "manual"
    if (is.null(names(dfs))) {
        ids <- paste0("species", seq_len(length(dfs)))
    } else {
        ids <- names(dfs)
    }
    stana@genes <- dfs
    stana@ids <- ids
    return(stana)
}


#' GF
#' @param dfs named list of dfs
#' @export
GF <- function(dfs) {
    stana <- new("stana")
    stana@type <- "manual"
    if (is.null(names(dfs))) {
        ids <- paste0("species", seq_len(length(dfs)))
    } else {
        ids <- names(dfs)
    }
    stana@kos <- dfs
    stana@ids <- ids
    return(stana)
}

#' fromHumann
#' import from HUMAnN profiles
#' @param df data frame (row.names as taxonomy and column names as sample names)
#' @export
#' @return stana object
fromHumann <- function(df) {
    ls <- row.names(df)
    tax <- ls[!startsWith(ls, "UNGROUPED")]
    tax <- tax[grepl("\\|", tax)]
    whole <- tax %>% strsplit("\\|") %>% vapply("[", 2, FUN.VALUE="a") %>% unique()
    whole <- whole[whole!="unclassified"]
    gfs <- lapply(whole, function(w) {
        tmp <- df[tax[endsWith(tax, w)], ]
        row.names(tmp) <- row.names(tmp) %>% strsplit("\\|") %>% sapply("[", 1)
        tmp
    }) %>% `setNames`(whole)
    stana <- new("stana")
    stana@kos <- gfs
    stana@ids <- names(gfs)
    stana@names <- names(gfs)
    stana@type <- "manual"
    stana
}


#' getGenes (concatenated to checkProfile)
#' 
#' Obtain gene matrix from midas merge directory from MIDAS.
#' Additionally, heatmap can be stored with clustering information.
#' 
#' @param midas_merge_dir output directory of merge_midas.py
#' @param candidate candidate species, default to all
#' @param pa "presabs" or "copynum"
#' @param km kmeans number if perform heatmap clustering
#' @param filUp filter parameter for upper limit
#' @param filDown filter parameter for lower limit
#' @param heatmap whether to draw heatmap
#' @param seed heatmap seed
#' @import GetoptLong
#' @importFrom methods new
#' @export
# getGenes <- function(midas_merge_dir,
#                      candidate="all",
#                      pa="presabs", km=20,
#                      filUp=1,filDown=0,
#                      heatmap=FALSE,
#                      seed=1) {
#     set.seed(seed)
#     mg <- new("stana")
#     dirLs <- list.files(midas_merge_dir)
#     specNames <- NULL
#     for (d in dirLs) {
#         if (dir.exists(paste0(midas_merge_dir,"/",d))){
#             specNames <- c(specNames, d)
#         }
#     }
#     # qqcat("@{specNames}\n")
#     mg@IDs <- specNames
#     mg@paFilterUp <- filUp
#     mg@paFilterDown <- filDown
#     if (candidate=="all") {
#       candSps <- specNames
#     } else {
#       candSps <- candidate
#     }
#     for (sp in candSps) {
#         qqcat("@{sp}\n")
#         df <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_",pa,".txt"), sep="\t",
#                          row.names=1, header=1)
#         qqcat("  Original gene count: @{dim(df)[1]}\n")
#         if (pa=="presabs") {
#           qqcat("  Gene present in all samples: @{sum(apply(df,1,sum)==dim(df)[2])}\n")
#         }
#         qqcat("  Gene absent in all samples: @{sum(apply(df,1,sum)==0)}\n")
#         if (pa=="presabs") {
#           presFil <- apply(df, 1, sum)>=dim(df)[2]*filUp
#           absFil <- apply(df, 1, sum)<=dim(df)[2]*filDown
#           filtDf <- df[!(presFil | absFil),]
#         } else {
#           qqcat("  For copy number, no filters were applied\n")
#           filtDf <- df
#         }
#         qqcat("  Filtered gene count: @{dim(filtDf)[1]}\n")
#         mg@Mat[[sp]] <- filtDf
#         if (heatmap){
#           qqcat("  Drawing heatmap\n")
#           hm <- draw(Heatmap(filtDf,
#                              show_row_names = FALSE,
#                              name=pa,
#                              row_km = km))
#           mg@Heatmap[[sp]] <- hm
#           dend <- row_dend(hm)
#           rowCl <- row_order(hm)
#           mg@geneCluster[[sp]] <- lapply(rowCl,
#                                      function(x) row.names(filtDf)[x])
#         }
#     }
#     mg
# }