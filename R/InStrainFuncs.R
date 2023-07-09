#' genomeHeatmap
#'
#' Draw a heatmap representing columns described in genome-wide
#' comparison table.
#' 
#' @param stana stana object
#' @param species interesting species
#' @param column column to plot, default to conANI
#' @param cl grouping, if NULL, it automatically obtain grouping from stana object
#' @return Heatmap by ComplexHeatmap
genomeHeatmap <- function(stana, species, column="conANI", cl=NULL) {
  if (stana@type!="InStrain") {stop("Please provide InStrain profiles")}
  if (is.null(cl)) {cl <- stana@cl}
  subMat <- a@genomeWideCompare[grepl(species,
                                      a@genomeWideCompare$genome),] |> 
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
  Heatmap(m, name=column,
          column_split=gr[colnames(m)])
}