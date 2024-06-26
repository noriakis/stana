#'
#' consensusSeq
#'
#' Output consensus sequences from merged SNV output.
#' Optionally, return phylogenetic tree inferred by `phangorn`.
#' If specified cluster of samples, additionally returns plot by `ggtree`.
#'
#' @param stana  stana object
#' @param species species vectors
#' @param argList parameters, passed to corresponding functions
#' @export
#' @rdname consensusseq
#'
consensusSeq <- function(stana,
	species=NULL, argList=list()){
    if (is.null(species)) {species <- stana@clearSnps}
    if (length(species)==0) {stop("No species available")}
	argList[["stana"]] <- stana
	argList[["species"]] <- species
	if (stana@type=="MIDAS2") {
		do.call("consensusSeqMIDAS2", argList)
	} else if (stana@type=="MIDAS1"){
		do.call("consensusSeqMIDAS1", argList)
	} else {
        cat("Calling general function\n")
		do.call("consensusSeqGeneral", argList)
	}
}

# #' consensusSeqMIDAS1
# #' 
# #' @param stana stana obj
# #' @param species candidate species vector
# #' @param mean_depth parameter for filtering
# #' @param fract_cov parameter for filtering
# #' @param site_depth parameter for filtering
# #' @param site_ratio parameter for filtering
# #' @param site_maf parameter for filtering
# #' @param allele_support parameter for filtering
# #' @param site_prev parameter for filtering
# #' @param cl cluster, if plot cladogram
# #' @param max_sites currently not implemented
# #' @param keep_samples currently not implemented
# #' @param exclude_samples currently not implemented
# #' @param rand_samples currently not implemented
# #' @param tree if perform tree inference using dist.ml()
# #' and NJ() in default parameters
# #' @param max_samples currently not implemented
# #' @importFrom phangorn read.phyDat dist.ml NJ
# #' @import ggtree ggplot2
# #' @importFrom phangorn read.phyDat
# #' @export
# consensusSeqMIDAS1Slow <- function(
# 	stana,
# 	species=NULL,
# 	mean_depth=0,
# 	fract_cov=0,
# 	site_depth=5,
# 	site_ratio=5.0,
# 	site_maf=0.01,
# 	allele_support=0.5,
# 	site_prev=0.9,
# 	cl=NULL,
# 	max_sites=Inf,
# 	tree=FALSE,
#     max_samples=Inf,
#     keep_samples=NULL,
#     exclude_samples=NULL,
#     rand_samples=NULL) {
#     ## site-list is currently not supported.
#     if (is.null(species)) {species <- stana@ids}
#     midas_merge_dir <- stana@mergeDir
# 	files <- c("depth","info","summary")
# 	retList <- list()
# 	if (sum(species %in% names(stana@snps))!=length(species)){
# 		stop("Species not included in loaded SNP tables")
# 	}
# 	for (sp in species) {
# 		qqcat("Beginning @{sp}\n")
# 		SPECIES <- list()
# 		if (is.null(stana@snpsInfo[[sp]])) {
# 			file <- "info"
# 			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_",file,".txt")
# 			SPECIES[[file]] <- read.table(filePath, sep="\t", header=1)			
# 		} else {
# 			SPECIES[["info"]] <- stana@snpsInfo[[sp]]
# 		}
# 		if (is.null(stana@snpsDepth[[sp]])) {
# 			file <- "depth"
# 			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_",file,".txt")
# 			SPECIES[[file]] <- read.table(filePath, sep="\t", header=1)	
# 		} else {
# 			SPECIES[["depth"]] <- stana@snpsDepth[[sp]]
# 		}


# 		if (dim(stana@snpsSummary)[1]==0) {
# 			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_summary.txt")
# 			snpsSummary <- read.table(filePath, header=1)
# 	        SPECIES[["summary"]] <- snpsSummary			
# 		} else {
# 			SPECIES[["summary"]] <- subset(stana@snpsSummary, stana@snpsSummary$species_id==sp)
# 		}			

# 		SPECIES[["freq"]] <- stana@snps[[sp]]
# 		siteNum <- dim(SPECIES[["freq"]])[1]
# 		qqcat("  Site number: @{siteNum}\n")
# 		SAMPLES <- list()

#         for (i in seq_len(nrow(SPECIES[["summary"]]))) {
#         	info <- SPECIES[["summary"]][i,]
#         	if (info$fraction_covered < fract_cov) {
#         		next
#         	} else if (info$mean_coverage < mean_depth) {
#         		next
#         	} else {
# 	        	SAMPLES[[info$sample_id]] <- list(mean_depth=info$mean_coverage,
# 	        		fract_cov=info$fraction_covered)
#         	}
#         }

#         for (sample in names(SAMPLES)) {
#         	SAMPLES[[sample]][["freq"]] <- as.numeric(SPECIES[["freq"]][sample][,1])
#         	SAMPLES[[sample]][["depth"]] <- as.numeric(SPECIES[["depth"]][sample][,1])
#         	## Append sample-wise filter
# 	        SAMPLES[[sample]][["filter"]][["site_depth"]] <- 
# 	            SAMPLES[[sample]][["depth"]] >=  site_depth
# 	        SAMPLES[[sample]][["filter"]][["depth_ratio"]] <- 
# 	            (SAMPLES[[sample]][["depth"]] / SAMPLES[[sample]]$mean_depth) <= site_ratio
# 	        SAMPLES[[sample]][["filter"]][["allele_support"]] <- 
# 	            apply(cbind(SAMPLES[[sample]][["freq"]], 1-SAMPLES[[sample]][["freq"]]), 1, max) >= allele_support
# 	    }
#         SITEFILTERS <- list()
#         retainedSites <- 0
#         for (i in seq_len(nrow(SPECIES[["freq"]]))) {
#         	if (retainedSites >= max_sites) {break}
#         	keepSamples <- NULL
#         	pooledMaf <- NULL
#         	for (sample in names(SAMPLES)) {
#         		keep <- sum(SAMPLES[[sample]][["filter"]][["site_depth"]][i],
#         		SAMPLES[[sample]][["filter"]][["depth_ratio"]][i],
#         		SAMPLES[[sample]][["filter"]][["allele_support"]][i])
#         		if (keep==3) {
#         			keepSamples <- c(keepSamples, sample)
#         			## Weight function is currently not supported.
#         			pooledMaf <- c(pooledMaf, SAMPLES[[sample]][["freq"]][i])
#         		}
#         	}
# 			if(is.null(keepSamples)) {pooledMaf <- 0}
#         	SITEFILTERS[[as.character(i)]][["prevalance"]] <- (length(keepSamples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev))
#         	SITEFILTERS[[as.character(i)]][["pooled_maf"]] <- mean(pooledMaf) >= site_maf
# 	        SITEFILTERS[[as.character(i)]][["ref_allele"]] <- SPECIES[["info"]][i,]$ref_allele %in% c("A","T","G","C")
# 	        SITEFILTERS[[as.character(i)]][["site_type"]] <- SPECIES[["info"]][i,]$site_type %in% c("1D","2D","3D","4D")
# 	        SITEFILTERS[[as.character(i)]][["locus_type"]] <- SPECIES[["info"]][i,]$locus_type %in% "CDS"

# 	        keepSite <- 
# 	        SITEFILTERS[[as.character(i)]][["prevalance"]]+
# 	        SITEFILTERS[[as.character(i)]][["pooled_maf"]]+
# 	        SITEFILTERS[[as.character(i)]][["ref_allele"]]+
# 	        SITEFILTERS[[as.character(i)]][["site_type"]]+
# 	        SITEFILTERS[[as.character(i)]][["locus_type"]]

# 	        if (keepSite==5) {
# 	        	retainedSites <- retainedSites + 1
# 	        	for (sample in names(SAMPLES)){
# 	        		if (!sample %in% keepSamples) {
# 	        			ap <- "-"
# 	        		} else if (SAMPLES[[sample]][["depth"]][i]==0) {
# 	        			ap <- "-"
# 	        		} else if (SAMPLES[[sample]][["freq"]][i]>=0.5) {
# 	        			ap <- SPECIES[["info"]][i,]$minor_allele
# 	        		} else {
# 	        			ap <- SPECIES[["info"]][i,]$major_allele
# 	        		}
# 	        		SAMPLES[[sample]][["consensus"]] <- 
# 	        		paste0(SAMPLES[[sample]][["consensus"]],ap)
# 	        	}
# 	        }
#         }
#         faName <- paste0(sp,"_consensus.fasta")
#         qqcat("  Outputting consensus sequence to @{faName}\n")
#         fileConn <- file(faName, open="w")
#         for (sample in names(SAMPLES)) {
# 			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
# 			cat(SAMPLES[[sample]][["consensus"]], file=fileConn, append=TRUE, sep="\n")
# 		}
# 		close(fileConn)
# 		tre <- read.phyDat(faName, format = "fasta")
# 		stana@fastaList[[sp]] <- tre
# 		if (tree) {
# 			dm <- dist.ml(tre, "F81")
# 			tre <- NJ(dm)
# 			stana@treeList[[sp]] <- tre
# 			if (!is.null(cl)) {
# 			    tre <- groupOTU(tre, cl)
# 			    tp <- ggtree(tre, aes(color=.data$group),
# 		               layout='circular',branch.length = "none") + # Return cladogram by default
# 			           geom_tippoint(size=3) + ggtitle(sp)+scale_color_manual(values=stana@colors)
# 			    stana@treePlotList[[sp]] <- tp
# 			}
# 		}
# 	}
# 	stana
# }


# #' consensusSeqMIDAS2
# #' @param stana stana obj
# #' @param species candidate species vector
# #' @param mean_depth parameter for filtering
# #' @param fract_cov parameter for filtering
# #' @param site_depth parameter for filtering
# #' @param site_ratio parameter for filtering
# #' @param site_maf parameter for filtering
# #' @param allele_support parameter for filtering
# #' @param site_prev parameter for filtering
# #' @param cl cluster, if plot cladogram
# #' @param max_sites default to Inf
# #' @param keep_samples currently not implemented
# #' @param exclude_samples currently not implemented
# #' @param rand_samples currently not implemented
# #' @param tree if perform tree inference
# #' @param max_samples currently not implemented
# #' @param verbose output current status
# #' @param output_seq whether to output actual FASTA file
# #' @export
# consensusSeqMIDAS2Slow <- function(
# 	stana,
# 	species,
# 	output_seq=FALSE,
# 	mean_depth=0,
# 	fract_cov=0,
# 	site_depth=5,
# 	site_ratio=5.0,
# 	site_maf=0.01,
# 	allele_support=0.5,
# 	site_prev=0.9,
# 	cl=NULL,
# 	max_sites=Inf,
# 	tree=FALSE,
#     max_samples=Inf,
#     keep_samples=NULL,
#     exclude_samples=NULL,
#     rand_samples=NULL,
#     verbose=FALSE) {
#     ## site-list is currently not supported.
# 	files <- c("depth","info")
# 	midas_merge_dir <- stana@mergeDir
# 	retList <- list()
# 	for (sp in species) {
# 		qqcat("Beginning calling for @{sp}\n")
# 		SPECIES <- list()
# 		## Load info and depth file per species
# 		if (is.null(stana@snpsInfo[[sp]])) {
# 			file <- "info"
# 			if (verbose) {
# 				qqcat("  Loading @{file]\n")
# 			}
# 			filePath <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_",file,".tsv.lz4")
#             filePathUn <- gsub(".lz4","",filePath)
# 	        system2("lz4", args=c("-d","-f",
# 	                                paste0(getwd(),"/",filePath),
# 	                                paste0(getwd(),"/",filePathUn)),
# 	                  stdout=FALSE, stderr=FALSE)
# 			SPECIES[[file]] <- read.table(filePathUn, row.names=1, header=1)
# 	        unlink(paste0(getwd(),"/",filePathUn))
# 		} else {
# 			SPECIES[["info"]] <- stana@snpsInfo[[sp]]
# 		}
# 		if (stana@type=="MIDAS2") {
# 			if (is.null(stana@snpsDepth[[sp]])) {
# 				file <- "depth"
# 				if (verbose) {
# 					qqcat("  Loading @{file]\n")
# 				}
# 				filePath <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_",file,".tsv.lz4")
# 	            filePathUn <- gsub(".lz4","",filePath)
# 		        system2("lz4", args=c("-d","-f",
# 		                                paste0(getwd(),"/",filePath),
# 		                                paste0(getwd(),"/",filePathUn)),
# 		                  stdout=FALSE, stderr=FALSE)
# 				SPECIES[[file]] <- read.table(filePathUn, row.names=1, header=1)
# 		        unlink(paste0(getwd(),"/",filePathUn))
# 			} else {
# 				SPECIES[["depth"]] <- stana@snpsDepth[[sp]]
# 			}			
# 		}
# 		SPECIES[["freqs"]] <- stana@snps[[sp]]
# 		siteNum <- dim(SPECIES[["freqs"]])[1]
# 		qqcat("  Site number: @{siteNum}\n")

# 		if (stana@type=="MIDAS2") {
# 			if (dim(stana@snpsSummary)[1]==0) {
# 				filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
# 				snpsSummary <- read.table(filePath, header=1)
# 		        SPECIES[["summary"]] <- subset(snpsSummary, snpsSummary$species_id==sp)			
# 			} else {
# 				SPECIES[["summary"]] <- subset(stana@snpsSummary, stana@snpsSummary$species_id==sp)
# 			}			
# 		}

# 		SAMPLES <- list()
#         for (i in seq_len(nrow(SPECIES[["summary"]]))) {
#         	if (verbose) {
#         		qqcat("  Filtering based on summary table using `fract_cov` and `mean_depth`\n")
#         	}
#         	info <- SPECIES[["summary"]][i,]
#         	if (verbose) {
# 	        	qqcat("    Processing @{info$sample_name}")        		
#         	}
#         	if (info$fraction_covered < fract_cov) {
#         		next
#         	} else if (info$mean_coverage < mean_depth) {
#         		next
#         	} else {
#         		if (verbose) {
# 	        		qqcat("  ... included, mean_coverage: @{info$mean_coverage}, fraction_covered: @{info$fraction_covered}\n")
#         		}
# 	        	SAMPLES[[info$sample_name]] <- list(mean_depth=info$mean_coverage,
# 	        		fract_cov=info$fraction_covered)
#         	}
#         }
#         for (sample in names(SAMPLES)) {
#         	if (verbose) {
# 	        	qqcat("    Processing @{sample} based on `site_depth`, `depth_ratio`, and `allel_support`\n")
#         	}
#         	SAMPLES[[sample]][["freqs"]] <- as.numeric(SPECIES[["freqs"]][sample][,1])
#         	SAMPLES[[sample]][["depth"]] <- as.numeric(SPECIES[["depth"]][sample][,1])
#         	## Append sample-wise filter
# 	        SAMPLES[[sample]][["filter"]][["site_depth"]] <- 
# 	            SAMPLES[[sample]][["depth"]] >=  site_depth
# 	        SAMPLES[[sample]][["filter"]][["depth_ratio"]] <- 
# 	            (SAMPLES[[sample]][["depth"]] / SAMPLES[[sample]]$mean_depth) <= site_ratio
# 	        SAMPLES[[sample]][["filter"]][["allele_support"]] <- 
# 	            apply(cbind(SAMPLES[[sample]][["freqs"]], 1-SAMPLES[[sample]][["freqs"]]), 1, max) >= allele_support
# 	    }
#         SITEFILTERS <- list()
#         retainedSites <- 0
#         for (i in seq_len(nrow(SPECIES[["freqs"]]))) {
#         	if (verbose) {qqcat("Site: @{i}\n")}
#         	if (retainedSites >= max_sites) {break}
#         	keepSamples <- NULL
#         	pooledMaf <- NULL
#         	for (sample in names(SAMPLES)) {
#         		if (verbose) {qqcat("  @{sample}")}
#         		keep <- sum(SAMPLES[[sample]][["filter"]][["site_depth"]][i],
#         		SAMPLES[[sample]][["filter"]][["depth_ratio"]][i],
#         		SAMPLES[[sample]][["filter"]][["allele_support"]][i])
#         		if (keep==3) {
#         			if (verbose) {
#         				qqcat(" ... Passed\n")
#         			}
#         			keepSamples <- c(keepSamples, sample)
#         			## Weight function is currently not supported.
#         			pooledMaf <- c(pooledMaf, SAMPLES[[sample]][["freqs"]][i])
#         		} else {
#         			if (verbose) {
#         				qqcat(" ... Not Passed\n")
#         			}
#         		}
#         	}
# 			if(is.null(keepSamples)) {pooledMaf <- 0}
#         	SITEFILTERS[[as.character(i)]][["prevalance"]] <- (length(keepSamples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev))
#         	SITEFILTERS[[as.character(i)]][["pooled_maf"]] <- mean(pooledMaf) >= site_maf
# 	        # SITEFILTERS[[as.character(i)]][["ref_allele"]] <- SPECIES[["info"]][i,]$ref_allele %in% c("A","T","G","C")
# 	        SITEFILTERS[[as.character(i)]][["site_type"]] <- SPECIES[["info"]][i,]$site_type %in% c("1D","2D","3D","4D")
# 	        SITEFILTERS[[as.character(i)]][["locus_type"]] <- SPECIES[["info"]][i,]$locus_type %in% "CDS"

# 	        keepSite <- 
# 	        SITEFILTERS[[as.character(i)]][["prevalance"]]+
# 	        SITEFILTERS[[as.character(i)]][["pooled_maf"]]+
# 	        SITEFILTERS[[as.character(i)]][["site_type"]]+
# 	        SITEFILTERS[[as.character(i)]][["locus_type"]]

# 	        if (keepSite==4) {
# 	        	retainedSites <- retainedSites + 1
# 	        	for (sample in names(SAMPLES)){
# 	        		if (!sample %in% keepSamples) {
# 	        			ap <- "-"
# 	        		} else if (SAMPLES[[sample]][["depth"]][i]==0) {
# 	        			ap <- "-"
# 	        		} else if (SAMPLES[[sample]][["freqs"]][i]>=0.5) {
# 	        			ap <- SPECIES[["info"]][i,]$minor_allele
# 	        		} else {
# 	        			ap <- SPECIES[["info"]][i,]$major_allele
# 	        		}
# 	        		SAMPLES[[sample]][["consensus"]] <- 
# 	        		paste0(SAMPLES[[sample]][["consensus"]],ap)
# 	        	}
# 	        }
#         }
#         faName <- paste0(sp,"_consensus.fasta")
#         qqcat("  Outputting consensus sequence to @{faName}\n")
        
#         fileConn<-file(faName, open="w")
#         for (sample in names(SAMPLES)) {
# 			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
# 			cat(SAMPLES[[sample]][["consensus"]], file=fileConn, append=TRUE, sep="\n")
# 		}
# 		close(fileConn)
# 		tre <- read.phyDat(faName, format = "fasta")
# 		stana@fastaList[[sp]] <- tre
# 		if (!output_seq) {
# 			unlink(faName)
# 		}
# 		if (tree) {
# 			dm <- dist.ml(tre, "F81")
# 			tre <- NJ(dm)
# 			stana@treeList[[sp]] <- tre
# 			if (!is.null(cl)) {
# 			    tre <- groupOTU(tre, cl)
# 			    tp <- ggtree(tre, aes(color=.data$group),
# 		               layout='circular',branch.length = "none") + # Return cladogram by default
# 			           geom_tippoint(size=3) + ggtitle(sp)+scale_color_manual(values=stana@colors)
# 			    stana@treePlotList[[sp]] <- tp
# 			}
# 		}
# 	}
# 	stana
# }