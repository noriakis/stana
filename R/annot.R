#' anno_PATRIC_keywords
#'
#' Annotate the heatmap using functional terms derived from PATRIC.
#' The function is derived and modified from the example of 
#' KeywordsEnrichment (https://github.com/jokergoo/KeywordsEnrichment) 
#' and using simplifyEnrichment package.
#' Can be replaced by providing simplyfyEnrichment::anno_word_cloud() 
#' with split and results from checkPATRICSimple().
#'
#' @param split named list of cluster 
#' @param genes gene names
#' @param fnc "pathway_name" or "ec_description"
#' @param removeHigh remove frequent words
#' @param removeAdditional remove these words
#'        passed to anno_word_cloud
#' @param argList will be passed to anno_word_cloud
#' @export
anno_PATRIC_keywords <- function(split, genes, fnc="pathway_name",
  removeHigh=TRUE, removeAdditional=NULL, argList=list()) {

  gene_list <- split(genes, split)
  gene_list <- lapply(gene_list, function(x) {
    x[!is.na(x)]
  })

  qqcat("Obtaining gene information from PATRIC server\n")
  lt <- checkPATRICSimple(gene_list, fnc)    
  names(lt) <- names(gene_list)

  if (removeHigh) {
    highTerms <- c("metabolism","biosynthesis","degradation")
  } else {
    highTerms <- NULL
  }
  argList[["align_to"]] <- split
  argList[["term"]] <- lt
  argList[["exclude_words"]] <- c(highTerms,removeAdditional)
  do.call("anno_word_cloud", argList)
}


#' anno_eggNOG_keywords
#'
#' Annotate the heatmap using functional terms derived from eggNOG.
#' The function is derived and modified from the example of 
#' KeywordsEnrichment (https://github.com/jokergoo/KeywordsEnrichment) 
#' and using simplifyEnrichment package.
#' Can be replaced by providing simplyfyEnrichment::anno_word_cloud()
#'
#' @param split named list of cluster 
#' @param genes gene names
#' @param tib result of checkEGGNOG()
#' @param removeHigh remove frequent words
#' @param removeAdditional remove these words
#'        passed to anno_word_cloud
#' @param argList will be passed to anno_word_cloud
#' @export
anno_eggNOG_keywords <- function(split, genes, tib,
  removeHigh=TRUE, removeAdditional=NULL, argList=list()) {

  gene_list <- split(genes, split)
  gene_list <- lapply(gene_list, function(x) {
    x[!is.na(x)]
  })

  lt <- lapply(gene_list, function(x) {
    (tib |> dplyr::filter(ID %in% x))$V2
  })

  names(lt) <- names(gene_list)

  if (removeHigh) {
    highTerms <- c("metabolism","biosynthesis","degradation","pathways","metabolic")
  } else {
    highTerms <- NULL
  }
  argList[["align_to"]] <- split
  argList[["term"]] <- lt
  argList[["exclude_words"]] <- c(highTerms,removeAdditional)
  do.call("anno_word_cloud", argList)
}