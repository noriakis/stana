#' exportInteractive
#' 
#' export the current stana object to interactive application 
#' for the convenient analysis and visualization
#' if `cl` and `meta` slot is filled, `meta` slot is exported.
#' 
#' @param stana stana object of type MIDAS2
#' @param out output directory
#' @param db db used to profile 'uhgg' or 'gtdb'
#' @param calcko calculate KO abundance
#' @param calctree calculate consensus tree
#' @param species candidate species, default to all the species
#' @param dataset_name dataset name
#' @export
#' @return output the results to specified directory
exportInteractive <- function(stana, out=".", db="uhgg", calcko=FALSE,
	calctree=FALSE, dataset_name=NULL,
	species=NULL) {
	if (is.null(dataset_name)) {
		dataset_name <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
	}
	dir.create(paste0(out,"/data"))
	dir.create(paste0(out,"/data/",dataset_name))

	# if (stana@type!="MIDAS2") {stop("This feature is for MIDAS2 only")}
	if (is.null(species)) {
		species <- stana@ids	
	}
	if (calcko) {
		for (sp in species) {
	        if (is.null(stana@kos[[sp]])) {
	        	if (!is.null(stana@eggNOG[[sp]])) {
		        	qqcat("Summarizing abundances @{sp}\n")
		            ko_tbl <- summariseAbundance(stana,sp = sp,
		                        checkEGGNOG(annot_file=stana@eggNOG[[sp]],
		                        	"KEGG_ko"),
		                        how="mean")
		            stana@kos[[sp]] <- ko_tbl	        		
	        	}
	        } else {
	            qqcat("Using pre-computed KO table\n")
	        }
	    }
	} else {
		for (sp in species) {
			if (is.null(stana@treeList[[sp]])) {
				if (is.null(stana@genes[[sp]])) {
					cat("Need at least the gene copy number table if tree is not available\n")
				}
			}
	        if (is.null(stana@kos[[sp]])) {
	        	## Insert gene table instead
	        	stana@kos[[sp]] <- stana@genes[[sp]]
	        }
	    }
	}
	
	for (tre in species) {
		if (is.null(stana@treeList[[tre]])) {
			qqcat("No tree for @{tre}\n")
			if (calctree) {
				qqcat("Calculating ... \n")
				stana <- consensusSeqMIDAS2(stana, tre)
                stana <- inferAndPlotTree(stana, tre)
			}
		}
	}
	if (length(stana@treeList)!=0) {
		treLen <- sum(sapply(stana@treeList, function(x) !is.null(x)))	
	} else {
		treLen <- 0
	}
	koLen <- length(stana@kos)

    all_samples_in_dataset1 <- unique(lapply(stana@treeList, function(x) {x$tip.label}) %>% unlist())
    all_samples_in_dataset2 <- unique(lapply(stana@kos, function(x) {colnames(x)}) %>% unlist())
    all_samples_in_dataset <- union(all_samples_in_dataset1, all_samples_in_dataset2)
	
	## If no metadata is available in meta slot:
    if (dim(stana@meta)[1]==0) {
        if (length(stana@cl)==0) {
            meta <- data.frame(all_samples_in_dataset) %>% `colnames<-`(c("samples"))
            meta[["label"]] <- meta[,1]
            row.names(meta) <- meta[,1]
        } else {
            meta <- NULL
            for (i in names(stana@cl)) {
                meta <- rbind(meta,
                    cbind(stana@cl[[i]],
                        rep(i, length(stana@cl[[i]]))))
            }
            meta <- meta |> data.frame() |> `colnames<-`(c("label","group"))
            row.names(meta) <- meta[["label"]]            
        }
    } else {
        meta <- stana@meta
    }
	
	qqcat("Tree number: @{treLen}, KO (or gene) number: @{koLen}\n")
	qqcat("Exporting ... \n")
	if (!dir.exists(out)) {
		dir.create(out)
	}
	dir.create(paste0(out,"/data/",dataset_name,"/tree"))
	for (i in names(stana@treeList)) {
		if (i %in% species) {
			tre <- stana@treeList[[i]]
			if (!is.null(tre)) {
				ape::write.tree(tre, paste0(out,"/data/",dataset_name,"/tree/",i,".cons.tree"))
			}			
		}
	}
	dir.create(paste0(out,"/data/",dataset_name,"/KOs"))
	for (i in names(stana@kos)) {
		if (i %in% species) {
			ko <- stana@kos[[i]]
			if (!is.null(ko)) {
				write.table(ko, paste0(out,"/data/",dataset_name,"/KOs/", i, ".txt"), sep="\t", quote=FALSE)						
			}
		}
	}

	write.table(meta, paste0(out,"/data/",dataset_name,"/metadata.tsv"), sep="\t", quote=FALSE)

    ## Copy the main script and run the app
    # file.copy(system.file("extdata", "app.R", package = "stana"), out)
    # shiny::runApp(paste0(out,"/app.R"))
	return(stana)
}