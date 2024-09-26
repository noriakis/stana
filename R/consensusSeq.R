#'
#' consensusSeq
#'
#' Output consensus sequences from merged SNV output.
#' Optionally, return phylogenetic tree inferred by `phangorn`.
#' If specified cluster of samples, additionally returns plot by `ggtree`.
#'
#' @param stana  stana object
#' @param species species vectors
#' @param argList parameters, passed to corresponding functions
#' @export
#' @rdname consensusseq
#' @return stana object
#'
consensusSeq <- function(stana,
	species=NULL, argList=list()){
    if (is.null(species)) {species <- stana@clearSnps}
    if (length(species)==0) {stop("No species available")}
	argList[["stana"]] <- stana
	argList[["species"]] <- species
	if (stana@type=="MIDAS2") {
		do.call("consensusSeqMIDAS2", argList)
	} else if (stana@type=="MIDAS1"){
		do.call("consensusSeqMIDAS1", argList)
	} else {
        cat("Calling general function\n")
		do.call("consensusSeqGeneral", argList)
	}
}
