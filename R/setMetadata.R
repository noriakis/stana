#' setMetadata
#' 
#' @param stana stana
#' @param meta metadata (data.frame)
#' row.name as ID of the trees
#' @export
#' @return stana object
#' 
setMetadata <- function(stana, meta) {
    meta[["id"]] <- row.names(meta)
    stana@meta <- meta
    return(stana)
}