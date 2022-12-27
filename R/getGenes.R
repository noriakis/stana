setClass("stana", slots=list(
                            type="character",
                            mergeDir="character",
                            ids="character",
                            snps="list",
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
                            heatmap="list",
                            geneCluster="list",
                            paFilterUp="numeric",
                            paFilterDown="numeric"))
setMethod("show",
  signature(object="stana"),
  function(object) {
    qqcat("Type: @{object@type}\n")
    qqcat("Species: @{length(object@ids)}\n")
    qqcat("Loaded SNV table: @{length(object@snps)}\n")
    qqcat("Loaded gene table (@{object@geneType}): @{length(object@genes)}\n")
  })
#' getGenes
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
#' @import ComplexHeatmap
#' @importFrom methods new
#' @importFrom ComplexHeatmap Heatmap
#' @export
getGenes <- function(midas_merge_dir,
                     candidate="all",
                     pa="presabs", km=20,
                     filUp=1,filDown=0,
                     heatmap=FALSE,
                     seed=1) {
    set.seed(seed)
    mg <- new("stana")
    dirLs <- list.files(midas_merge_dir)
    specNames <- NULL
    for (d in dirLs) {
        if (dir.exists(paste0(midas_merge_dir,"/",d))){
            specNames <- c(specNames, d)
        }
    }
    # qqcat("@{specNames}\n")
    mg@IDs <- specNames
    mg@paFilterUp <- filUp
    mg@paFilterDown <- filDown
    if (candidate=="all") {
      candSps <- specNames
    } else {
      candSps <- candidate
    }
    for (sp in candSps) {
        qqcat("@{sp}\n")
        df <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_",pa,".txt"), sep="\t",
                         row.names=1, header=1)
        qqcat("  Original gene count: @{dim(df)[1]}\n")
        if (pa=="presabs") {
          qqcat("  Gene present in all samples: @{sum(apply(df,1,sum)==dim(df)[2])}\n")
        }
        qqcat("  Gene absent in all samples: @{sum(apply(df,1,sum)==0)}\n")
        if (pa=="presabs") {
          presFil <- apply(df, 1, sum)>=dim(df)[2]*filUp
          absFil <- apply(df, 1, sum)<=dim(df)[2]*filDown
          filtDf <- df[!(presFil | absFil),]
        } else {
          qqcat("  For copy number, no filters were applied\n")
          filtDf <- df
        }
        qqcat("  Filtered gene count: @{dim(filtDf)[1]}\n")
        mg@Mat[[sp]] <- filtDf
        if (heatmap){
          qqcat("  Drawing heatmap\n")
          hm <- draw(Heatmap(filtDf,
                             show_row_names = FALSE,
                             name=pa,
                             row_km = km))
          mg@Heatmap[[sp]] <- hm
          dend <- row_dend(hm)
          rowCl <- row_order(hm)
          mg@geneCluster[[sp]] <- lapply(rowCl,
                                     function(x) row.names(filtDf)[x])
        }
    }
    mg
}