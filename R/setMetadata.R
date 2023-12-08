#' setMetadata
#' 
#' @param stana stana
#' @param meta metadata (data.frame)
#' row.name as ID of the trees
#' 
#' @return stana object
#' 
setMetadata(stana, meta) {
    meta[["id"]] <- row.names(meta)
    stana@meta <- meta
    return(stana)
}