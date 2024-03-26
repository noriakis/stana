#' plotKEGGPathway
#'
#' plot the KEGG pathway
#' The function uses `cl` slot and output statistics based on the groups.
#' If the groups number is above two, the sum value across the group 
#' will be used (not using grouping information).
#'
#' @param stana stana class object
#' @param species candidate species
#' @param pathway_id KEGG PATHWAY id
#' @param multi_scale if TRUE, use ``ggh4x`` to plot multiple scales per species
#' @param how how to concatenate abundance per KO
#' @param eggnog use eggNOG annotation in stana object
#' @param eps pseudovalue added when claculating log2 fold change
#' @param kegg_name_match passed to node_numeric function
#' @param color_list color list if multi_scale=TRUE,
#' passed to scale_fill_multi()
#' @param only_ko only calculates KO and return stana
#' @param cl grouping (named list). If not specified, use stana@cl slot
#' @param summarize summarize the value if multiple species are set ("+")
#' @param point_mode default to FALSE, show whether the corresponding KOs
#' are present in species by points
#' @param sp_colors species colors used in point_mode
#' @param nudge_y nudge the y positions in point_mode
#' @param point_size point size in point_mode
#' @export
#' @return list of plots or plot
#' @importFrom ggkegg pathway
#' @importFrom ggplotify as.ggplot
#'
plotKEGGPathway <- function(stana, species, pathway_id,
                            cl=NULL, multi_scale=FALSE,
                            eggnog=TRUE, kegg_name_match="all",
                            how=mean, eps=1e-2, color_list=NULL,
                            only_ko=FALSE, summarize=FALSE, nudge_y=0,
                            point_mode=FALSE, sp_colors=NULL,
                            point_size=1, statMethod="t"){
    sum_flag <- FALSE	
    if (is.null(sp_colors)) {
    	## Assign random colors
    	sp_colors <- scales::brewer_pal(palette="RdBu")(length(species)) |>
		setNames(species)
	}
    if (length(species)==1 & summarize) {stop("summarize option is intended for multiple species")}
    if (is.null(color_list)) {
        color_list <- list(
            scales::brewer_pal(palette = "YlGnBu")(6),
            scales::brewer_pal(palette = "RdPu")(6),
            scales::brewer_pal(palette = "YlOrRd")(6)
        )
    }
    if (length(color_list[[1]])<length(species)) stop("not sufficient color list")
    if (sum(species %in% names(stana@genes))!=length(species)) {stop("Not all the species are available in genes slot")}
    if (is.null(cl)) cl <- stana@cl
    if (length(names(cl))>2) {GetoptLong::qqcat("Group number above two, changing to sum value\n");
        sum_flag <- TRUE}
    if (eggnog) {
        if (sum(species %in% names(stana@eggNOG))!=length(species)) stop("Not all the annotation for species are available")
    }
  
    spnum <- length(species)
  
    lfcs <- list()
    for (sp in species) {
        if (is.null(stana@kos[[sp]])) {
            ko_tbl <- summariseAbundance(stana,sp = sp,
                        checkEGGNOG(annot_file=stana@eggNOG[[sp]], "KEGG_ko"),
                        how=how)
            stana@kos[[sp]] <- ko_tbl      
        } else {
            qqcat("Using pre-computed KO table\n")
        }
    }

    ## Return stana object only
    if (only_ko) return(stana)


    if (summarize) {
        commons <- Reduce(intersect, lapply(stana@kos[species], function(x) row.names(x)))
        commoncol <- Reduce(intersect, lapply(stana@kos[species], function(x) colnames(x)))
        ko_tbl <- Reduce("+", lapply(stana@kos[species], function(x) x[commons,commoncol]))

        if (sum_flag) {
          lfcs[["Sum"]] <- apply(ko_tbl, 1, sum)
        } else {
            qqcat("@{sp}: @{names(cl)[1]} / @{names(cl)[2]}\n")
            lfcs[[sp]] <- L2FC(ko_tbl, cl[[1]], cl[[2]], method=statMethod)
        }
    } else {
        for (sp in species) {
            ko_tbl <- stana@kos[[sp]]
            if (sum_flag) {
                lfcs[[sp]] <- apply(ko_tbl, 1, sum)
            } else {
                qqcat("@{sp}: @{names(cl)[1]} / @{names(cl)[2]}\n")
                lfcs[[sp]] <- L2FC(ko_tbl, cl[[1]], cl[[2]], method=statMethod)
            }
        }        
    }


    ## Obtain graph
    
    
    
    if (point_mode) {
    	qqcat("Point mode enabled\n")
        graphList <- list()
        for (pid in pathway_id) {
            g <- ggkegg::pathway(pid)
            for (sp in species) {
	            g <- g |>
	            	dplyr::mutate(!!sp := ggkegg::node_numeric(lfcs[[sp]]))        	
            }
            nds <- g |> tidygraph::activate("nodes") |> data.frame(check.names=FALSE)            
            number <- rowSums(!is.na(nds[,species]))
            g <- g |> dplyr::mutate(number=number) |> dplyr::mutate(space=width/(number-1))
            nds <- g |> tidygraph::activate("nodes") |> data.frame(check.names=FALSE)            
            
            spnum <- length(species)

			gg <- ggraph(g, layout="manual", x=x, y=y)
			
			for (sp in seq_along(species)) {
    			## Drop the NA values
                nds_tmp <- nds[!is.na(nds[[species[sp]]]),]
                num <- sp-1
                nds_tmp <- nds_tmp |> dplyr::filter(type=="ortholog") |> 
                    dplyr::mutate(tmp_x=xmin+ifelse(number==1, width/2, space*!!num))
    			gg <- gg + geom_node_point(
    				        aes(fill="transparent"),
                            x=nds_tmp$tmp_x,
                            y=nds_tmp$y+nudge_y,
                            color=sp_colors[species[sp]],
                            size=point_size,
                            data=nds_tmp)
    		}
    		plotter <- gg+ggkegg::overlay_raw_map()+theme_void()
	
			## Make pseudo-legend
			pseudo <- data.frame(x=rep(1, spnum),
			           y=rep(1, spnum),
			           group=species)
			dd2 <- ggplot(pseudo, aes(x=x, y=y))+
			  geom_node_point(aes(color=group))+
			  scale_color_manual(values=sp_colors, name="Species")+
			  guides(color = guide_legend(override.aes = list(shape=19, size=5)))+
			  theme_void()
			grobs <- ggplot_gtable(ggplot_build(dd2))
			num <- which(sapply(grobs$grobs, function(x) x$name) == "guide-box")
			legendGrob <- grobs$grobs[[num]]
			overlaidGtable <- ggplot_gtable(ggplot_build(plotter))
			num2 <- which(sapply(overlaidGtable$grobs, function(x) x$name) == "guide-box")
			overlaidGtable$grobs[[num2]] <- legendGrob
			returng <- ggplotify::as.ggplot(overlaidGtable)
            graphList[[pid]] <- returng
		}
        return(graphList)
    }
    

    if (summarize) {
        graphList <- list()
        for (pid in pathway_id) {
            g <- ggkegg::pathway(pid)
            g <- g |> dplyr::mutate(Sum := ggkegg::node_numeric(lfcs[[sp]], name="name",how=kegg_name_match))
        }
        V(g)$space <- V(g)$width
        nds <- g |> tidygraph::activate("nodes") |> data.frame()
        nds <- nds[nds$type %in% "ortholog",]
        gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)
      
        gg <- gg + ggkegg::geom_node_rect(
            aes(fill=Sum,
                filter=.data$type %in% "ortholog")
            )
        gg <- gg + ggkegg::overlay_raw_map()+ggplot2::theme_void()
        graphList[[pid]] <- gg
    } else {
        graphList <- list()
        for (pid in pathway_id) {
          g <- ggkegg::pathway(pid)
          for (sp in species) {
            g <- g |> 
              dplyr::mutate(!!sp := ggkegg::node_numeric(lfcs[[sp]], name="name",how=kegg_name_match))
          }
          V(g)$space <- V(g)$width/spnum
          nds <- g |> tidygraph::activate("nodes") |> data.frame()
          nds <- nds[nds$type %in% "ortholog",]
          gg <- ggraph(g, layout="manual", x=.data$x, y=.data$y)
          if (!multi_scale) {
            for (i in seq_len(spnum)) {
              nudge <- i-1
              gg <- gg + ggkegg::geom_node_rect(
                aes(fill=!!sym(species[i]),
                    filter=.data$type %in% "ortholog"),
                xmin=nds$xmin+nds$space*nudge,
                xmax=nds$xmin+i*nds$space
              )
            }
            gg <- gg + ggkegg::overlay_raw_map()+ggplot2::theme_void()
            graphList[[pid]] <- gg
          } else {
            ## Add new scale recursively ##
            for (i in seq_len(spnum)) {
              nudge <- i-1
              gg <- gg + ggkegg::geom_node_rect(
                aes(fill= !!sym(species[i]),
                    filter=.data$type %in% "ortholog"),
                xmin=nds$xmin+nds$space*nudge,
                xmax=nds$xmin+i*nds$space,
              )
              gg$layers[[i]]$mapping[[species[i]]] <- gg$layers[[i]]$mapping[["fill"]]
              gg$layers[[i]]$mapping[["fill"]] <- NULL
            }
            gg <- gg + ggkegg::overlay_raw_map()+ggplot2::theme_void()
            sps <- list()
            for (i in seq_len(length(species))) {
              sps[[i]] <- species[i]
            }
            gg <- gg +
              ggh4x::scale_fill_multi(
                aesthetics = species,
                name = sps,
                colours = color_list[seq_len(length(species))],
                guide = guide_colorbar(barheight = unit(50, "pt"))
              )
          }
          graphList[[pid]] <- gg
        }        
    }
    if (length(graphList)==1) {
      return(graphList[[pid]])
    } else {
      return(graphList)
    }
}