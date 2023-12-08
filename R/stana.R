#' Note that gene summary and depth is not loaded in MIDAS and MIDAS2.
#' @import methods
#' @importFrom utils object.size
#' @importFrom grDevices colorRampPalette
#' @slot snpsInfo snp info
#' @slot snpsDepth snp depth
#' @slot snpsSummary snp summary
#' @slot kos stores ko table
#' @slot type type of pipeline, such as `MIDAS`
#' @slot eggNOG list of path for eggNOG mapper v2 results for species
setClass("stana", slots=list(
                            type="character",
                            cl="list",
                            meta="data.frame",
                            mergeDir="character",
                            db="character",
                            ids="character",
                            comparisonTable="data.table",
                            genomeWideCompare="data.table",
                            strainClusters="data.table",
                            snps="list",
                            eggNOG="list",
                            kos="list",
                            snpsInfo="list",
                            snpsDepth="list",
                            snpsSummary="data.frame",
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
                            colors="vector",
                            geneCluster="list",
                            paFilterUp="numeric",
                            paFilterDown="numeric"))
setMethod("show",
  signature(object="stana"),
  function(object) {
    qqcat("Type: @{object@type}\n")
    qqcat("Directory: @{object@mergeDir}\n")
    qqcat("Species number: @{length(object@ids)}\n")
    if (object@type %in% c("MIDAS","MIDAS2")) {
      qqcat("Filter type: @{object@sampleFilter}, number: @{object@sampleFilterVal}, proportion: @{object@sampleFilterPer}\n")
    }
    if (length(object@snps)!=0) {
      qqcat("Loaded SNV table: @{length(object@snps)}\n")      
    }
    if (object@type %in% c("MIDAS","MIDAS2")) {
      qqcat("  Species cleared SNV filter: @{length(object@clearSnps)}\n")
    }
    if (length(object@genes)!=0) {
      qqcat("Loaded gene table (@{object@geneType}): @{length(object@genes)}\n")
    }
    if (object@type %in% c("MIDAS","MIDAS2")) {
      qqcat("  Species cleared gene filter: @{length(object@clearGenes)}\n")
    }
    print(object.size(object), units="auto")
  })



#' @export
setGeneric("getSlot",
    function(x, ...) standardGeneric("getSlot"))

setMethod("getSlot", "stana",
    function(x, slot) attr(x, slot))

#' @export
setGeneric("getFasta",
    function(x) standardGeneric("getFasta"))

setMethod("getFasta", "stana",
    function(x) attr(x, "fastaList"))

#' @export
setGeneric("getTree",
    function(x) standardGeneric("getTree"))

setMethod("getTree", "stana",
    function(x) attr(x, "treeList"))

#' @export
setGeneric("getID",
    function(x) standardGeneric("getID"))

setMethod("getID", "stana",
    function(x) attr(x, "ids"))

#' @export
setGeneric("getAdonis",
    function(x) standardGeneric("getAdonis"))

setMethod("getAdonis", "stana",
    function(x) attr(x, "adonisList"))

#' @export
setGeneric("getCl",
    function(x) standardGeneric("getCl"))

setMethod("getCl", "stana",
    function(x) attr(x, "cl"))

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