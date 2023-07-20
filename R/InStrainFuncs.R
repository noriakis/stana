#' genomeHeatmap
#'
#' Draw a heatmap representing columns described in genome-wide
#' comparison table.
#' 
#' @param stana stana object
#' @param species interesting species
#' @param column column to plot, default to conANI
#' @param cl grouping, if NULL, it automatically obtain grouping from stana object
#' @param heatmapArgs list of parameters passed to Heatmap()
#' @return Heatmap by ComplexHeatmap
#' @export
genomeHeatmap <- function(stana, species, column="conANI", cl=NULL, heatmapArgs=list()) {
  if (stana@type!="InStrain") {stop("Please provide InStrain profiles")}
  if (is.null(cl)) {cl <- stana@cl}
  subMat <- stana@genomeWideCompare[grepl(species,
                                      stana@genomeWideCompare$genome),] |> 
    data.frame()
  subMat <- subMat[,c("name1","name2",column)]
  edge_mat <- subMat |>
    tidyr::pivot_wider(names_from=name2,values_from=!!column) |>
    data.frame(check.names=FALSE)
  gr <- NULL
  if (length(cl)>0) {
    for (cln in names(cl)) {
      gr <- c(gr, rep(cln, length(cl[[cln]])) |> setNames(cl[[cln]]))
    }
  }
  edge_mat$name1 <- NULL
  m <- edge_mat
  m[lower.tri(m)] <- t(m)[lower.tri(t(m))]
  heatmapArgs[["matrix"]] <- m
  heatmapArgs[["name"]] <- column
  heatmapArgs[["column_split"]] <- gr[colnames(m)]
  do.call(Heatmap, heatmapArgs)
}



#' strainClusterHeatmap
#'
#' Draw a heatmap representing presence/absence of strain cluster
#' 
#' 
#' @param stana stana object
#' @param species interesting species
#' @param cl grouping, if NULL, it automatically obtain grouping from stana object
#' @param heatmapArgs list of parameters passed to Heatmap()
#' @return Heatmap by ComplexHeatmap
#' @export
strainClusterHeatmap <- function(stana, species, cl=NULL, heatmapArgs=list()) {
  if (stana@type!="InStrain") {stop("Please provide InStrain profiles")}
  if (is.null(cl)) {cl <- stana@cl}
  sc <- stana@strainClusters
  candsc <- sc[grepl(species, sc$genome),]

  subMat <- candsc[,c("cluster","sample")]
  subMat$present <- 1
  edge_mat <- subMat |>
      tidyr::pivot_wider(names_from=sample,values_from=present) |>
      data.frame(check.names=FALSE)
  gr <- NULL
  if (length(cl)>0) {
      for (cln in names(cl)) {
          gr <- c(gr, rep(cln, length(cl[[cln]])) |> setNames(cl[[cln]]))
      }
  }
  row.names(edge_mat) <- edge_mat$cluster
  edge_mat$cluster <- NULL
  edge_mat[is.na(edge_mat)] <- 0
  heatmapArgs[["matrix"]] <- edge_mat
  heatmapArgs[["name"]] <- "cluster_present"
  heatmapArgs[["column_split"]] <- gr[colnames(edge_mat)]
  do.call(Heatmap, heatmapArgs)
}