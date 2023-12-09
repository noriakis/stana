#' exportInteractive
#' 
#' export the current stana object to interactive application 
#' for the convenient analysis and visualization
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
exportInteractive <- function(stana, out=".", db="uhgg", calcko=TRUE,
	calctree=TRUE, dataset_name=NULL,
	species=NULL) {
	if (is.null(dataset_name)) {
		dataset_name <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
	}
	dir.create(paste0(out,"/data"))
	dir.create(paste0(out,"/data/",dataset_name))

	if (stana@type!="MIDAS2") {stop("This feature is for MIDAS2 only")}
	if (is.null(species)) {
		species <- stana@ids	
	}
	if (calcko) {
		for (sp in species) {
	        if (is.null(stana@kos[[sp]])) {
	        	qqcat("Summarizing abundances @{sp}\n")
	            ko_tbl <- summariseAbundance(stana,sp = sp,
	                        checkEGGNOG(annot_file=stana@eggNOG[[sp]],
	                        	"KEGG_ko"),
	                        how="mean")
	            stana@kos[[sp]] <- ko_tbl      
	        } else {
	            qqcat("Using pre-computed KO table\n")
	        }
	    }
	}
	
	for (tre in species) {
		if (is.null(stana@treeList[[tre]])) {
			qqcat("No tree for @{tre}\n")
			if (calctree) {
				qqcat("Calculating ... \n")
				stana <- consensusSeqMIDAS2(stana, tre, tree=TRUE,
					max_sites=50) # temp	
			}
		}
	}
	treLen <- sum(sapply(stana@treeList, function(x) !is.null(x)))
	koLen <- length(stana@kos)
	meta <- NULL
	for (i in names(stana@cl)) {
		meta <- rbind(meta,
			cbind(stana@cl[[i]],
				rep(i, length(stana@cl[[i]]))))
	}
	meta <- meta |> data.frame() |> `colnames<-`(c("label","group"))
	qqcat("Tree number: @{treLen}, KO number: @{koLen}\n")
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
			write.table(ko, ppaste0(out,"/data/",dataset_name,"/KOs/", i, ".txt"), sep="\t", quote=FALSE)			
		}
	}
	write.table(meta, paste0(out,"/data/",dataset_name,"/metadata.tsv"), sep="\t", quote=FALSE)
	return(stana)
}