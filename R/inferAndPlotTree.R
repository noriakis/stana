#' inferAndPlotTree
#' 
#' infer the tree and plot.
#' If there is already a tree in the list from `treeList` slot, inferring will not be performed.
#' The function uses dist.ml and NJ for inference from FASTA format, and otherwise uses dist and NJ
#' for the inference. if metadata is available, use them to plot metadata around the tree.
#' 
#' @param stana stana object
#' @param species species to plot
#' If NULL, first species in fasta list is assigned
#' @param cl optional, cluster to plot
#' @param model dist.ml model
#' @param tree_args passed to dist function
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
#' @param deleteZeroDepth delete zero depth position
#' @param treeFun tree inferring function in phangorn
#' if `FastTree` is specified, FastTree will be called.
#' In this case, FastTree should be in PATH.
#' @param tree_only return tree only instead of stana object
#' @param subset_samples subset samples (for matrix only)
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggstar geom_star
#' @import phangorn
#' @importFrom scico scale_fill_scico_d scale_fill_scico scale_color_scico scale_color_scico_d
#' 
#' @export
inferAndPlotTree <- function(stana, species=NULL, cl=NULL,
	dist_method="dist.ml", meta=NULL, layout="circular",
	target="fasta", IDs=NULL, use_point=FALSE, branch_col="black",
	tree_args=list(), branch.length="none", point_size=2,
    subset_samples=NULL,
    deleteZeroDepth=TRUE, treeFun="upgma", tree_only=FALSE) {
	if (is.null(cl)) {cl <- stana@cl}
	if (!is.null(meta)) {
		meta <- checkMeta(stana, meta)
	}
	if (is.null(species)) {species <- names(stana@fastaList)}
	if (!(target %in% c("fasta","snp","gene","KO"))) {stop("Please specify appropriate target (snp, gene, fasta)")}
	if (target!="fasta" & treeFun=="FastTree") {
		stop("FastTree should be used with FASTA target")
	}
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
            if (!is.null(subset_samples)) {
                mat <- mat[, intersect(colnames(mat),subset_samples)]
            }
			tree_args[["x"]] <- t(mat)
			dm <- do.call(dist, tree_args)
		} else if (target=="KO") {
            mat <- stana@kos[[sp]]
            if (!is.null(IDs)) {
                mat <- mat[intersect(row.names(mat),IDs), ]
            }
            if (!is.null(subset_samples)) {
                mat <- mat[, intersect(colnames(mat),subset_samples)]
            }
            tree_args[["x"]] <- t(mat)
            dm <- do.call(dist, tree_args)      
        } else {
			mat <- stana@snps[[sp]]
            if (deleteZeroDepth) {
                mat <- mat[rowSums(mat==-1)==0, ]
                cat("Position number:", dim(mat)[1], "\n")      
            } else {
            	mat[ mat == -1] <- NA
            }
			if (!is.null(IDs)) {
				mat <- mat[intersect(row.names(mat),IDs), ]
			}
            if (!is.null(subset_samples)) {
                mat <- mat[, intersect(colnames(mat),subset_samples)]
            }
			tree_args[["x"]] <- t(mat)
			dm <- do.call(dist, tree_args)			
		}
        ## Difficult to call do.call
        ## If needed, FastTree can be called.
        ## In this case, FastTree should be in PATH
        if (treeFun=="FastTree") {
        	## Write out the FASTA to the current dir
        	fa <- stana@fastaList[[sp]]
        	nam <- paste0(sp, "_consensus_MSA_stana_tmp.fa")
        	if (file.exists(nam)) {
        		cat("File already exists!\n")
        	} else {
	        	write.phyDat(fa, nam, format="fasta")    		
        	}
        	trenam <- paste0(sp, "_consensus_tree_stana_tmp.tree")
        	if (file.exists(trenam)) {
        		cat("Tree file already exists!\n")
        	}
            system2("FastTree", args=c("-out", trenam, "-nt", nam),
                stdout=TRUE, stderr=TRUE)
            tre <- read.tree(trenam)
        	## Call FastTree
        	## Read tree and save
        }
        if (treeFun=="upgma") {
            tre <- upgma(dm)   
        } else {
            tre <- NJ(dm)
        }
        if (tree_only) {
            return(tre)
        }
        ## Negative edge length would be present in NJ
        ## Reference: https://boopsboops.blogspot.com/2010/10/negative-branch-lengths-in-neighbour.html
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
            
            
            if (branch_col %in% colnames(stana@meta)) {
                ## NA value is fixed
                tmp_class <- class(stana@meta[[branch_col]])
                stana@meta$label <- stana@meta$id
                addNow <- stana@meta[,c(branch_col, "label")]
                tre <- full_join(tre, addNow, by="label")
                p <- ggtree(tre, layout=flyt, branch.length=branch.length, mapping=aes(color=.data[[branch_col]]))
                if (tmp_class %in% c("integer","numeric","logical")) { ## T/F is treated as numeric
                    p <- p + scico::scale_color_scico(na.value="grey", palette=sample(scico::scico_palette_names(), 1))
                } else {
                    p <- p + scico::scale_color_scico_d(na.value="grey", palette=sample(scico::scico_palette_names(), 1))
                }
            } else {
                p <- ggtree(tre, layout=flyt, branch.length=branch.length, color=branch_col)            
            }
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
	                        mapping = aes(y=id, fill=.data[[tmp_show_cv]])
	                    )                	
	                } else {
	                    p <- p + geom_fruit(
	                        data = meta,
	                        geom = geom_star,
	                        size = point_size,
	                        mapping = aes(y=id, fill=.data[[tmp_show_cv]]),
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