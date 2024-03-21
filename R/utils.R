#' Discretize copy number matrix at specified cutoff
#' @param stana stana object
#' @param species species id
#' @param cutoff cutoff value, default to 0.35
#' @export
cnDiscretize <- function(stana, species, cutoff=0.35) {
	df <- stana@genes[[species]]
	apply(df, 2, function(x) ifelse(x>cutoff, 1, 0)) |> data.frame()
}

#' setTree
#' @param stana stana object
#' @param species species ID
#' @param tre tree to be set
#' @export
#' @return stana
setTree <- function(stana, species, tre) {
	if (class(tre)!="phylo") {stop("Please provide phylo object")}
	stana@treeList[[species]] <- tre
	return(stana)
}


#' setGroup
#' @param stana stana object
#' @param cl grouping list
#' @export
#' @return stana
setGroup <- function(stana, cl) {
	stana@cl <- cl
	return(stana)
}


#' returnGenes
#' return corresponding genes queried by  SNV position
#' @param stana stana object
#' @param species species name
#' @param snvs snv name
#' @export
#' @return named vector
returnGenes <- function(stana, species, snvs) {
  check <- stana@snpsInfo[[species]][snvs, ]
  check <- check[check$gene_id!="None",]
  if (dim(check)[1]>0) {
    return(check$gene_id |> setNames(row.names(check)))
  } else {
    stop("No genes mapped")
  }
}

#' @noRd
listToNV <- function(l) {
    nm <- NULL
    val <- NULL
    for (i in names(l)) {
        val <- c(val, l[[i]])
        nm <- c(nm, rep(i, length(l[[i]])))
    }
    names(nm) <- val
    return(nm)
}


#' plotSNVinfo
#' 
#' @param stana stana object
#' @param sp candidate species
#' @export
#' @return ggplot
plotSNVInfo <- function(stana, sp) {
    d <- stana@snpsInfo[[sp]]
    tbl <- data.frame(table(paste0(d$major_allele,"/",d$minor_allele)))
    if ("locus_type" %in% colnames(d)) {
        cdsc <- data.frame(table(d$locus_type))
        cdsc <- paste(paste0(cdsc$Var1, " ", cdsc$Freq), collapse=" ")
        title <- paste0("Total: ",dim(d)," (", cdsc, ")")
    }
    ggplot(tbl, aes(x=Var1, y=Freq)) +
        geom_col()+
        xlab("major_allele/minor_allele")+
        cowplot::theme_cowplot()+
        ggtitle(title)
}

#' plotSNVSummary
#' @param stana stana object
#' @param sp species_id
#' @param param parameter to plot
#' @export 
#' @return ggplot object
plotSNVSummary <- function(stana, sp, param="mean_coverage") {
    df <- stana@snpsSummary
    if (dim(df)[1]==0) {stop("SNV summary not available")}
    df <- df %>% dplyr::filter(species_id == sp)
    if (length(stana@cl)!=0) {
        df[["group"]] <- listToNV(stana@cl)[df$sample_name]
        ggplot(df, aes(x=group, y=.data[[param]])) +
            geom_boxplot() + cowplot::theme_cowplot()
    } else {
        ggplot(df, aes(y=.data[[param]])) +
            geom_boxplot() + cowplot::theme_cowplot()
    }
}