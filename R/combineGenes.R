#' combineGenes
#' 
#' Use common genes across multiple datasets to
#' merge gene copy numbers
#' 
#' @param stana_list stana list
#' @param species species ID
#' @export
#' @return new stana object
combineGenes <- function(stana_list, species) {
  if (!is.list(stana_list)) {stop("Please provide list of stana object")}
  each_ID <- lapply(stana_list, function(x) x@genes[[species]] |> row.names())
  intersected <- Reduce(intersect, each_ID)
  if (length(intersected)==0) {stop("No common genes")}
  
  qqcat("Common genes: @{length(intersected)}\n")
  
  merged <- do.call(cbind, lapply(stana_list, function(x) x@genes[[species]][intersected,]))
  
  ## Metadata
  new_stana <- new("stana")
  new_stana@genes[[species]] <- merged
  cls <- lapply(stana_list, function(x) x@cl)
  cls <- do.call(c, cls)
  
  ## Warning if same label
  ovlgr <- length(Reduce(intersect, lapply(stana_list, function(x) names(x@cl))))
  if (ovlgr>0) {qqcat("Duplicate label found in group\n")}

  new_stana@ids <- species
  new_stana@type <- stana@type
  new_stana@cl <- cls
  new_stana@colors <- getColors(cls)

  return(new_stana)
}

