#' Discretize copy number matrix at specified cutoff
#' @param stana stana object
#' @param species species id
#' @param cutoff cutoff value, default to 0.35
#' @export
cnDiscretize <- function(stana, species, cutoff=0.35) {
	df <- stana@genes[["100224"]]
	apply(df, 2, function(x) ifelse(x>cutoff, 1, 0)) |> data.frame()
}