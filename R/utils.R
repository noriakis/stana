#' Discretize copy number matrix at specified cutoff
#' @param stana stana object
#' @param species species id
#' @param cutoff cutoff value, default to 0.35
#' @return discretized data.frame
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
	if (!is(tre, "phylo")) {stop("Please provide phylo object")}
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


#' scaler.NNLM
#' @param nmf NNLM object
#' @param target coef or basis
#' @export
#' @return list
scaler.NNLM <- function(nmf, target="coef") {
    if (target=="basis") {
        scaledW <- apply(nmf$W, 2, function(x) x / sum(x))
        Sii <- apply(nmf$W, 2, function(x) sum(x))
        scaledH <- apply(nmf$H, 2, function(x) {
            x * Sii
        })
        list("W"=scaledW, "H"=scaledH)        
    } else {
        scaledH <- apply(nmf$H, 1, function(x) x / sum(x)) %>% t()
        Sii <- apply(nmf$H, 1, function(x) sum(x))
        scaledW <- apply(nmf$W, 1, function(x) {
            x * Sii
        }) %>% t()
        list("W"=scaledW, "H"=scaledH)
    }

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
	if (!is.null(stana@includeSNVID[[sp]])) {cat_subtle("# The set SNV ID information is not used in this function\n")}
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
        scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
        ggtitle(title)
}

#' plotSNVSummary
#' @param stana stana object
#' @param sp species_id
#' @param param parameter to plot
#' @param perSample plot using no sample information
#' @export 
#' @return ggplot object
plotSNVSummary <- function(stana, sp, param="mean_coverage", perSample=FALSE) {
    df <- stana@snpsSummary
    if (dim(df)[1]==0) {stop("SNV summary not available")}
    df <- df %>% dplyr::filter(species_id == sp)
    if (perSample) {
        med <- df[[param]] %>% median(na.rm=TRUE)
        if (length(stana@cl)!=0) {
            df[["group"]] <- listToNV(stana@cl)[df$sample_name]
            ggplot(df, aes(x=sample_name, y=.data[[param]])) +
                geom_col(aes(fill=group)) +
                scale_fill_manual(values=stana@colors)+
                scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
                cowplot::theme_cowplot()+
                theme(axis.text.x = element_blank())+
                geom_hline(yintercept = med, lty=2)
        } else {
            ggplot(df, aes(x=sample_name, y=.data[[param]])) +
                geom_col() +
                scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
                cowplot::theme_cowplot()+
                theme(axis.text.x = element_blank())+
                geom_hline(yintercept = med, lty=2)
        }           
    } else {
        if (length(stana@cl)!=0) {
            df[["group"]] <- listToNV(stana@cl)[df$sample_name]
            ggplot(df, aes(x=group, y=.data[[param]])) +
                geom_boxplot(aes(fill=.data[[param]]), alpha=0.5) +
                scale_fill_manual(values=stana@colors)+
                cowplot::theme_cowplot()
        } else {
            ggplot(df, aes(y=.data[[param]])) +
                geom_boxplot() + cowplot::theme_cowplot()
        }        
    }

}