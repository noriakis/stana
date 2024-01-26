#' plotTree
#' 
#' infer the tree and plot
#' The function uses dist.ml and NJ for inference
#' if metadata is available, use them to plot metadata around the tree.
#' 
#' @param stana stana object
#' @param species species to plot
#' If NULL, first species in fasta list is assigned
#' @param cl optional, cluster to plot
#' @param model dist.ml model
#' @param tree_args passed to dist function in phangorn
#' @param dist_method dist method in phangorn, default to dist.ml
#' if target is not fasta, ordinally `dist` method wil be used.
#' @param branch.length branch length, default to "none", cladogram
#' @param meta default to NULL, column name in `meta` slot
#' @param point_size point size
#' @param layout layout to be used in ggtree
#' @param target fasta, snp, or gene. default to fasta.
#' @param IDs if IDs is provided and target is snps or genes,
#' subset the matrix to this ID.
#' @param use_point use geom_point for plotting discrete metadata,
#' instead of geom_star.
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggstar geom_star
#' @importFrom scico scale_fill_scico_d scale_fill_scico
#' 
#' @export
plotTree <- function(stana, species=NULL, cl=NULL,
	dist_method="dist.ml", meta=NULL, layout="circular",
	target="fasta", IDs=NULL, use_point=FALSE,
	tree_args=list(), branch.length="none", point_size=2) {
	if (is.null(cl)) {cl <- stana@cl}
	if (!is.null(meta)) {
		meta <- checkMeta(stana, meta)
	}
	if (is.null(species)) {species <- names(stana@fastaList)}
	if (!(target %in% c("fasta","snp","gene"))) {stop("Please specify appropriate target")}
	for (sp in species) {
		if (target=="fasta") {
			tre <- stana@fastaList[[sp]]
			
			## Infer tree
			tree_args[["x"]] <- tre
			dm <- do.call(dist_method, tree_args)
		} else if (target=="gene") {
			mat <- stana@genes[[sp]]
			if (!is.null(IDs)) {
				mat <- mat[intersect(row.names(mat),IDs), ]
			}
			tree_args[["x"]] <- t(mat)
			dm <- do.call(dist, tree_args)
		} else {
			mat <- stana@snps[[sp]]
			if (!is.null(IDs)) {
				mat <- mat[intersect(row.names(mat),IDs), ]
			}
			tree_args[["x"]] <- t(mat)
			dm <- do.call(dist, tree_args)			
		}
		tre <- NJ(dm)
		stana@treeList[[sp]] <- tre			
		
		## Plot tree		
		if (!is.null(meta)) {
            ## Add NA row to unknown label in tree
            ## not in metadata
            for (i in tre$tip.label) {
                if (!i %in% row.names(meta)) {
                    cur <- row.names(meta)
                    meta <- rbind(meta, rep(NA, ncol(meta)))
                    row.names(meta) <- c(cur, i)
                }
            }
            
            all_cols <- colnames(meta)
            show_meta <- all_cols[all_cols!="id"]
            
            for (tmp_show_cv in show_meta) {
            	tmp_class <- class(meta[[tmp_show_cv]])
                change_cv <- meta[[tmp_show_cv]]
                change_cv[change_cv==""] <- NA
                change_cv[change_cv=="-"] <- NA
                change_cv[is.infinite(change_cv)] <- NA

                if (tmp_class %in% c("integer","numeric")) {
                    meta[[tmp_show_cv]] <- as.numeric(change_cv)
                } else {
                    meta[[tmp_show_cv]] <- as.factor(change_cv)
                }
            }

            meta <- meta[tre$tip.label,]


            if (layout=="rectangular") {
                qqcat("Force setting layout to 'circular'\n")
                flyt <- "circular"
            } else {
                flyt <- layout
            }

            ## Select random palette from scico
            scp <- sample(scico::scico_palette_names(),
            	length(colnames(meta)), replace=FALSE)
            names(scp) <- colnames(meta)

            ## Select random shape from ggstar
            starsh <- sample(1:30, length(colnames(meta)), replace=FALSE)
            names(starsh) <- colnames(meta)

            p <- ggtree(tre, layout=flyt, branch.length=branch.length)
            for (tmp_show_cv in show_meta) {
                if (is.numeric(meta[[tmp_show_cv]])) {
                    ## If is numeric
                    p <- p + geom_fruit(
                        data = meta,
                        geom = geom_col,
                        mapping = aes(y=id, fill=.data[[tmp_show_cv]],
                        	x=.data[[tmp_show_cv]]),
                    ) +
                    scale_fill_scico(palette = scp[tmp_show_cv]) +
                    new_scale_fill()
                } else {
                	if (use_point) {
	                    p <- p + geom_fruit(
	                        data = meta,
	                        geom = geom_point,
	                        size = point_size, shape=21, 
	                        mapping = aes(y=id, x=0.5, fill=.data[[tmp_show_cv]])
	                    )                	
	                } else {
	                    p <- p + geom_fruit(
	                        data = meta,
	                        geom = geom_star,
	                        size = point_size,
	                        mapping = aes(y=id, x=0.5, fill=.data[[tmp_show_cv]]),
	                        starshape = starsh[tmp_show_cv]
	                    )
	                }
                    ## If is discrete
					p <- p + 
                    scale_fill_scico_d(palette = scp[tmp_show_cv], begin=0, end=0.5) +
                    new_scale_fill()
                }
            }
		    stana@treePlotList[[sp]] <- p
		} else {
		    tre <- groupOTU(tre, cl)
		    tp <- ggtree(tre,
		        layout=layout,branch.length = branch.length) + # Return cladogram by default
		        geom_tippoint(size=point_size, aes(color=.data$group)) +
		        ggtitle(sp)+
		        scale_color_manual(values=stana@colors)
		    stana@treePlotList[[sp]] <- tp			
		}
	}
    stana
}


#' @noRd
checkMeta <- function(stana, meta) {
	meta_col <- stana@meta |> colnames()
	int_cols <- intersect(meta_col, meta)
	int_len <- int_cols |> length()
	if (length(int_len)==0) {
		stop("No column name is available in the `meta` slot")
	}
	return(stana@meta[, c("id", int_cols)])
}