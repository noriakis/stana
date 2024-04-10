
#' consensusSeqGeneral
#' @param stana stana obj
#' @param species candidate species vector
#' @param cl cluster, if plot cladogram
#' @param keep_samples the samples to keep
#' @param tree if TRUE, perform tree inference
#' @param verbose output current status
#' @param output_seq whether to output actual FASTA file
#' @param return_mat return matrix of characters
#' @param allele_columns columns specifying major and minor allele
#' @export
consensusSeqGeneral <- function(
	stana,
	species=NULL,
	cl=NULL,
	keep_samples=NULL,
	output_seq=FALSE,
	tree=FALSE,
    return_mat=FALSE,
    allele_columns=c("major_allele","minor_allele"),
    verbose=FALSE) {
    
    if (is.null(species)) {species <- stana@clearSnps}
    if (length(species)==0) {stop("No species available")}
	retList <- list()
	for (sp in species) {
		cat_subtle("# Beginning calling for ", sp, "\n", sep="")
		SPECIES <- list()
		SPECIES[["freqs"]] <- stana@snps[[sp]]
		
		if (!is.null(stana@includeSNVID[[sp]])) {
			cat_subtle("# The set SNV ID information (",
				length(stana@includeSNVID[[sp]]), ") is used.\n")
			SPECIES[["freqs"]] <- SPECIES[["freqs"]][stana@includeSNVID[[sp]], ]
		}
		
		siteNum <- dim(SPECIES[["freqs"]])[1]
		qqcat("  Site number: @{siteNum}\n")
		
		
		if (is.null(stana@snpsInfo[[sp]])) {
			stop("Please provide data frame for snpsInfo slot")
		} else {
			if (length(intersect(allele_columns, colnames(stana@snpsInfo[[sp]])))<2) {
				stop("Snps info data frame has to have major_allele and minor_allele column")
			}
    		SPECIES[["info"]] <- data.frame(stana@snpsInfo[[sp]])
			if (!is.null(stana@includeSNVID[[sp]])) {
				SPECIES[["info"]] <- SPECIES[["info"]][stana@includeSNVID[[sp]], ]
			}    		
		}

        sp_minor_allele <- SPECIES[["info"]][[allele_columns[2]]]
        sp_major_allele <- SPECIES[["info"]][[allele_columns[1]]]
        
        allele_list <- lapply(seq_len(nrow(SPECIES[["freqs"]])), function(i) {
	        	if (!is.null(keep_samples)) {
	        		samples <- keep_samples
	        	} else {
	        		samples <- colnames(SPECIES[["freqs"]])
	        	}
        	    
		    	per_sample_allele <- lapply(samples, function(sample) { ## BOTTLENECK [3]
		    		if (SPECIES[["freqs"]][[sample]][i] == -1) {
		    		    return("-")
		    		}
                    if (SPECIES[["freqs"]][[sample]][i] >= 0.5) {
		    			return(sp_minor_allele[i])
		    		} else {
		    			return(sp_major_allele[i])
		    		}
		    	})
	    	unlist(per_sample_allele)        	
            })
        
        if (return_mat) {
        	mat <- do.call(cbind, allele_list)
        	row.names(mat) <- names(site_filters)
	        return(mat) 
        }
        
        allele_list <- apply(do.call(cbind, allele_list), 1, function(x) paste0(x, collapse="")) |>
        setNames(colnames(SPECIES[["freqs"]]))
        faName <- paste0(sp,"_consensus.fasta")
        if (output_seq) {
            qqcat("  Outputting consensus sequence to @{faName}\n")    	
        }
        
        fileConn<-file(faName, open="w")
        for (sample in colnames(SPECIES[["freqs"]])) {
			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
			cat(allele_list[sample], file=fileConn, append=TRUE, sep="\n")
		}
		close(fileConn)
		tre <- read.phyDat(faName, format = "fasta")
		stana@fastaList[[sp]] <- tre
		if (!output_seq) {
			unlink(faName)
		}
		if (tree) {
			dm <- dist.ml(tre, "F81")
			tre <- NJ(dm)
			stana@treeList[[sp]] <- tre
			if (!is.null(cl)) {
			    tre <- groupOTU(tre, cl)
			    tp <- ggtree(tre, aes(color=.data$group),
		               layout='circular',branch.length = "none") + # Return cladogram by default
			           geom_tippoint(size=3) + ggtitle(sp)+scale_color_manual(values=stana@colors)
			    stana@treePlotList[[sp]] <- tp
			}
		}
	}
	stana
}


#' consensusSeqMIDAS2Fast
#' @param stana stana obj
#' @param species candidate species vector
#' @param mean_depth parameter for filtering
#' @param fract_cov parameter for filtering
#' @param site_depth parameter for filtering
#' @param site_ratio parameter for filtering
#' @param site_maf parameter for filtering
#' @param allele_support parameter for filtering
#' @param site_prev parameter for filtering
#' @param cl cluster, if plot cladogram
#' @param max_sites default to Inf
#' @param keep_samples currently not implemented
#' @param exclude_samples currently not implemented
#' @param rand_samples currently not implemented
#' @param tree if perform tree inference
#' @param max_samples currently not implemented
#' @param verbose output current status
#' @param output_seq whether to output actual FASTA file
#' @param locus_type locus type to be included, default to CDS
#' @param site_list site list to be included
#' @param return_mat return matrix of characters
#' @export
consensusSeqMIDAS2 <- function(
	stana,
	species=NULL,
	output_seq=FALSE,
	mean_depth=0,
	fract_cov=0,
	site_depth=5,
	site_ratio=5.0,
	site_maf=0.01,
	allele_support=0.5,
	site_prev=0.9,
	locus_type="CDS",
	site_list=NULL,
	cl=NULL,
	max_sites=Inf,
	tree=FALSE,
    max_samples=Inf,
    keep_samples=NULL,
    exclude_samples=NULL,
    rand_samples=NULL,
    return_mat=FALSE,
    verbose=FALSE) {
    
    if (is.null(species)) {species <- stana@clearSnps}
    if (length(species)==0) {stop("No species available")}
	files <- c("depth","info")
	midas_merge_dir <- stana@mergeDir
	retList <- list()
	for (sp in species) {
		cat_subtle("# Beginning calling for ", sp, "\n", sep="")
		SPECIES <- list()
		## Load info and depth file per species
		
		if (!is.null(stana@includeSNVID[[sp]])) {
			cat_subtle("# The set SNV ID information (",
				length(stana@includeSNVID[[sp]]), ") is used.\n")
			site_list <- stana@includeSNVID[[sp]]
		}
		
		if (is.null(stana@snpsInfo[[sp]])) {
			file <- "info"
			if (verbose) {
				qqcat("  Loading @{file]\n")
			}
			filePath <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_",file,".tsv.lz4")
            filePathUn <- gsub(".lz4","",filePath)
	        system2("lz4", args=c("-d","-f",
	                                paste0(getwd(),"/",filePath),
	                                paste0(getwd(),"/",filePathUn)),
	                  stdout=FALSE, stderr=FALSE)
			SPECIES[[file]] <- read.table(filePathUn, row.names=1, header=1)
	        unlink(paste0(getwd(),"/",filePathUn))
		} else {
			SPECIES[["info"]] <- stana@snpsInfo[[sp]]
		}
		if (stana@type=="MIDAS2") {
			if (is.null(stana@snpsDepth[[sp]])) {
				file <- "depth"
				if (verbose) {
					qqcat("  Loading @{file]\n")
				}
				filePath <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_",file,".tsv.lz4")
	            filePathUn <- gsub(".lz4","",filePath)
		        system2("lz4", args=c("-d","-f",
		                                paste0(getwd(),"/",filePath),
		                                paste0(getwd(),"/",filePathUn)),
		                  stdout=FALSE, stderr=FALSE)
				SPECIES[[file]] <- read.table(filePathUn, row.names=1, header=1)
		        unlink(paste0(getwd(),"/",filePathUn))
			} else {
				SPECIES[["depth"]] <- stana@snpsDepth[[sp]]
			}			
		}
		SPECIES[["freqs"]] <- stana@snps[[sp]]
		siteNum <- dim(SPECIES[["freqs"]])[1]
		cat_subtle("# Original Site number: ", siteNum, "\n", sep="")

		if (stana@type=="MIDAS2") {
			if (dim(stana@snpsSummary)[1]==0) {
				filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
				snpsSummary <- read.table(filePath, header=1)
		        SPECIES[["summary"]] <- subset(snpsSummary, snpsSummary$species_id==sp)			
			} else {
				SPECIES[["summary"]] <- subset(stana@snpsSummary, stana@snpsSummary$species_id==sp)
			}			
		}

        SAMPLES <- lapply(SPECIES[["summary"]]$sample_name, function(x) {
        	info <- SPECIES[["summary"]][SPECIES[["summary"]]$sample_name == x, ]
        	if (!is.null(keep_samples)) {
            	if (!(x %in% keep_samples)) {
            		return(NULL)
            	}        		
        	}
        	if (info$fraction_covered < fract_cov) {
        		return(NULL)
        	} else if (info$mean_coverage < mean_depth) {
        		return(NULL)
        	} else {
        		if (verbose) {
	        		qqcat("  @{x} ... included, mean_coverage: @{info$mean_coverage}, fraction_covered: @{info$fraction_covered}\n")
        		}
	        	return(list(mean_depth=info$mean_coverage,
	        		fract_cov=info$fraction_covered))
        	}
        })
        cat_subtle("#  Profiled samples: ", length(SAMPLES), "\n", sep="")
        names(SAMPLES) <- SPECIES[["summary"]]$sample_name
        SAMPLES <- SAMPLES[vapply(SAMPLES, function(x) !is.null(x), FUN.VALUE=TRUE)]
        cat_subtle("#  Included samples: ", length(SAMPLES), "\n", sep="")
        
        if (verbose) {
            qqcat("Beginning site-wise filtering\n")
        }
        site_filters <- lapply(names(SAMPLES), function(sample) {
        	freqs <- SPECIES[["freqs"]][sample][,1] |> as.numeric()
        	depths <- SPECIES[["depth"]][sample][,1] |> as.numeric()
        	site_depth_filter <- depths >= site_depth
        	depth_ratio_filter <- (depths / SAMPLES[[sample]]$mean_depth) <= site_ratio
        	## freqs == -1 in allele_support        	        	
        	allele_support <- apply(cbind(freqs, 1-freqs), 1, max) >= allele_support ## BOTTLENECK [1]
        	allele_support[which(freqs==-1)] <- FALSE
        	summed <- site_depth_filter + depth_ratio_filter + allele_support
        	list(freqs, depths, site_depth_filter, depth_ratio_filter, allele_support, summed) |>
        	setNames(c("freqs","depths","site_depth","depth_ratio","allele_support","flag"))
        })
        names(site_filters) <- names(SAMPLES)

        sp_site_type <- SPECIES[["info"]]$site_type
        sp_locus_type <- SPECIES[["info"]]$locus_type
        sp_minor_allele <- SPECIES[["info"]]$minor_allele
        sp_major_allele <- SPECIES[["info"]]$major_allele
        
        keep_site_filter_list <- lapply(seq_len(nrow(SPECIES[["freqs"]])), function(i) {
        	keep_samples <- lapply(names(site_filters), function(sample) { ## BOTTLENECK [2]
        	    if (site_filters[[sample]][["flag"]][i]==3) {
        	    	list(sample, site_filters[[sample]][["freqs"]][i])
        	    } else {
        	    	NULL
        	    }
        	})
        	
        	keep_samples <- keep_samples[vapply(keep_samples, function(x) !is.null(x), FUN.VALUE=TRUE)]     	
        	kept_sample <- lapply(keep_samples, function(x) x[[1]]) |> unlist()
        	kept_sample_num <- length(kept_sample)
        	pooled_maf <- lapply(keep_samples, function(x) x[[2]]) |> unlist() |> mean()

			if(length(keep_samples)==0) {pooled_maf <- 0}
			
        	keep_site_filter <- 
        	    sum(c(
        	    	(length(keep_samples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev)),
        	        pooled_maf >= site_maf,
	                sp_site_type[i] %in% c("1D","2D","3D","4D"),
	                sp_locus_type[i] %in% locus_type))
	        return(keep_site_filter)
        }) |> unlist()
        
                
        if (!is.null(site_list)) {
        	cat_subtle("# site_list specified: ", length(site_list), "\n", sep="")
        	positions <- which(row.names(SPECIES[["info"]]) %in% site_list)
        	## Ignoring max_sites if site_list is provided.
        	## All the site in site_list is included, and filtered positions are 
        	## denoted as "-"
            allele_list <- lapply(positions, function(i) {
		    	per_sample_allele <- lapply(names(site_filters), function(sample) {
		    		if (site_filters[[sample]][["flag"]][i]!=3) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["depths"]][i] == 0) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["freqs"]][i] == -1) {
		    		    ap <- "-"	
		    		} else if (keep_site_filter_list[i]!=4) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["freqs"]][i] >= 0.5) {
		    			ap <- sp_minor_allele[i]
                    } else {
		    			ap <- sp_major_allele[i]
		    		}
		    		return(ap)
		    	})
	    	unlist(per_sample_allele)        	
            })
        } else {
	        positions <- which(keep_site_filter_list==4)
	        if (!is.infinite(max_sites)) {positions <- positions[1:max_sites]}
            allele_list <- lapply(positions, function(i) {
		    	per_sample_allele <- lapply(names(site_filters), function(sample) { ## BOTTLENECK [3]
		    		if (site_filters[[sample]][["flag"]][i]!=3) {
		    			return("-")
                    }
		    		if (site_filters[[sample]][["depths"]][i] == 0) {
		    			return("-")
                    }
		    		if (site_filters[[sample]][["freqs"]][i] == -1) {
		    		    return("-")
		    		}
                    if (site_filters[[sample]][["freqs"]][i] >= 0.5) {
		    			return(sp_minor_allele[i])
		    		} else {
		    			return(sp_major_allele[i])
		    		}
		    	})
	    	unlist(per_sample_allele)        	
            })
        }
        
        if (return_mat) {
        	mat <- do.call(cbind, allele_list)
        	row.names(mat) <- names(site_filters)
        	colnames(mat) <- row.names(SPECIES[["info"]])[positions]
	        return(mat) 
        }
        
        allele_list <- apply(do.call(cbind, allele_list), 1, function(x) paste0(x, collapse="")) |>
        setNames(names(site_filters))
        faName <- paste0(sp,"_consensus.fasta")
        if (output_seq) {
            qqcat("  Outputting consensus sequence to @{faName}\n")    	
        }
        
        fileConn<-file(faName, open="w")
        for (sample in names(site_filters)) {
			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
			cat(allele_list[sample], file=fileConn, append=TRUE, sep="\n")
		}
		close(fileConn)
		tre <- read.phyDat(faName, format = "fasta")
		stana@fastaList[[sp]] <- tre
		if (!output_seq) {
			unlink(faName)
		}
		if (tree) {
			dm <- dist.ml(tre, "F81")
			tre <- NJ(dm)
			stana@treeList[[sp]] <- tre
			if (!is.null(cl)) {
			    tre <- groupOTU(tre, cl)
			    tp <- ggtree(tre, aes(color=.data$group),
		               layout='circular',branch.length = "none") + # Return cladogram by default
			           geom_tippoint(size=3) + ggtitle(sp)+scale_color_manual(values=stana@colors)
			    stana@treePlotList[[sp]] <- tp
			}
		}
	}
	stana
}



#' consensusSeqMIDAS1
#' 
#' @param stana stana obj
#' @param species candidate species vector
#' @param mean_depth parameter for filtering
#' @param fract_cov parameter for filtering
#' @param site_depth parameter for filtering
#' @param site_ratio parameter for filtering
#' @param site_maf parameter for filtering
#' @param allele_support parameter for filtering
#' @param site_prev parameter for filtering
#' @param cl cluster, if plot cladogram
#' @param max_sites currently not implemented
#' @param keep_samples currently not implemented
#' @param exclude_samples currently not implemented
#' @param rand_samples currently not implemented
#' @param tree if perform tree inference using dist.ml()
#' and NJ() in default parameters
#' @param max_samples currently not implemented
#' @param verbose print output
#' @param output_seq output the FASTA file
#' @param return_mat return character matrix
#' @param locus_type locus type to be included
#' @param site_list site list to be included
#' @importFrom phangorn read.phyDat dist.ml NJ
#' @import ggtree ggplot2
#' @importFrom phangorn read.phyDat
#' @export
consensusSeqMIDAS1 <- function(
	stana,
	species=NULL,
	mean_depth=0,
	fract_cov=0,
	site_depth=5,
	site_ratio=5.0,
	site_maf=0.01,
	allele_support=0.5,
	site_prev=0.9,
	cl=NULL,
	max_sites=Inf,
	tree=FALSE,
	locus_type="CDS",
    max_samples=Inf,
    keep_samples=NULL,
    exclude_samples=NULL,
    rand_samples=NULL,
    verbose=FALSE,
    site_list=NULL,
    return_mat=FALSE,
    output_seq=FALSE) {
    ## site-list is currently not supported.
    if (is.null(species)) {species <- stana@clearSnps}
    if (length(species)==0) {stop("No species available")}
    midas_merge_dir <- stana@mergeDir
	files <- c("depth","info","summary")
	retList <- list()
	if (sum(species %in% names(stana@snps))!=length(species)){
		stop("Species not included in loaded SNP tables")
	}
	for (sp in species) {
		qqcat("Beginning @{sp}\n")
		SPECIES <- list()
		if (is.null(stana@snpsInfo[[sp]])) {
			file <- "info"
			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_",file,".txt")
			SPECIES[[file]] <- read.table(filePath, sep="\t", header=1)			
		} else {
			SPECIES[["info"]] <- stana@snpsInfo[[sp]]
		}
		if (is.null(stana@snpsDepth[[sp]])) {
			file <- "depth"
			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_",file,".txt")
			SPECIES[[file]] <- read.table(filePath, sep="\t", header=1)	
		} else {
			SPECIES[["depth"]] <- stana@snpsDepth[[sp]]
		}


		if (dim(stana@snpsSummary)[1]==0) {
			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_summary.txt")
			snpsSummary <- read.table(filePath, header=1)
	        SPECIES[["summary"]] <- snpsSummary			
		} else {
			SPECIES[["summary"]] <- subset(stana@snpsSummary, stana@snpsSummary$species_id==sp)
		}			

		SPECIES[["freqs"]] <- stana@snps[[sp]]
		siteNum <- dim(SPECIES[["freqs"]])[1]
		qqcat("  Site number: @{siteNum}\n")
		
        SAMPLES <- lapply(SPECIES[["summary"]]$sample_id, function(x) {
        	info <- SPECIES[["summary"]][SPECIES[["summary"]]$sample_id == x, ]
        	if (!is.null(keep_samples)) {
            	if (!(x %in% keep_samples)) {
            		return(NULL)
            	}        		
        	}

        	if (info$fraction_covered < fract_cov) {
        		return(NULL)
        	} else if (info$mean_coverage < mean_depth) {
        		return(NULL)
        	} else {
        		if (verbose) {
	        		qqcat("  ... included, mean_coverage: @{info$mean_coverage}, fraction_covered: @{info$fraction_covered}\n")
        		}
	        	return(list(mean_depth=info$mean_coverage,
	        		fract_cov=info$fraction_covered))
        	}
        })
        
        qqcat("  Profiled samples: @{length(SAMPLES)}\n")
        names(SAMPLES) <- SPECIES[["summary"]]$sample_id
        SAMPLES <- SAMPLES[vapply(SAMPLES, function(x) !is.null(x), FUN.VALUE=TRUE)]
        qqcat("  Included samples: @{length(SAMPLES)}\n")

        site_filters <- lapply(names(SAMPLES), function(sample) {
        	freqs <- SPECIES[["freqs"]][sample][,1] |> as.numeric()
        	depths <- SPECIES[["depth"]][sample][,1] |> as.numeric()
        	site_depth_filter <- depths >= site_depth
        	depth_ratio_filter <- (depths / SAMPLES[[sample]]$mean_depth) <= site_ratio
        	## freqs == -1 in allele_support        	        	
        	allele_support <- apply(cbind(freqs, 1-freqs), 1, max) >= allele_support
        	allele_support[which(freqs==-1)] <- FALSE
        	summed <- site_depth_filter + depth_ratio_filter + allele_support
        	list(freqs, depths, site_depth_filter, depth_ratio_filter, allele_support, summed) |>
        	setNames(c("freqs","depths","site_depth","depth_ratio","allele_support","flag"))
        })
        names(site_filters) <- names(SAMPLES)
        
        sp_site_type <- SPECIES[["info"]]$site_type
        sp_locus_type <- SPECIES[["info"]]$locus_type
        sp_minor_allele <- SPECIES[["info"]]$minor_allele
        sp_major_allele <- SPECIES[["info"]]$major_allele
        sp_ref_allele <- SPECIES[["info"]]$ref_allele
        
        keep_site_filter_list <- lapply(seq_len(nrow(SPECIES[["freqs"]])), function(i) {
        	keep_samples <- lapply(names(site_filters), function(sample) {
        	    if (site_filters[[sample]][["flag"]][i]==3) {
        	    	list(sample, site_filters[[sample]][["freqs"]][i])
        	    } else {
        	    	NULL
        	    }
        	})
        	
        	keep_samples <- keep_samples[vapply(keep_samples, function(x) !is.null(x), FUN.VALUE=TRUE)]     	
        	kept_sample <- lapply(keep_samples, function(x) x[[1]]) |> unlist()
        	kept_sample_num <- length(kept_sample)
        	pooled_maf <- lapply(keep_samples, function(x) x[[2]]) |> unlist() |> mean()

			if(length(keep_samples)==0) {pooled_maf <- 0}
			
        	keep_site_filter <- 
        	    sum(c(
        	    	(length(keep_samples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev)),
        	        pooled_maf >= site_maf,
        	        sp_ref_allele[i] %in% c("A","T","G","C"),
	                sp_site_type[i] %in% c("1D","2D","3D","4D"),
	                sp_locus_type[i] %in% locus_type))
	        return(keep_site_filter)
        }) |> unlist()
        
        ## 5 is good in MIDAS1
        if (!is.null(site_list)) {
        	positions <- which(row.names(SPECIES[["info"]]) %in% site_list)
        	## Ignoring max_sites if site_list is provided.
        	## All the site in site_list is included, and filtered positions are 
        	## denoted as "-"
            allele_list <- lapply(positions, function(i) {
		    	per_sample_allele <- lapply(names(site_filters), function(sample) {
		    		if (site_filters[[sample]][["flag"]][i]!=3) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["depths"]][i] == 0) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["freqs"]][i] == -1) {
		    		    ap <- "-"	
		    		} else if (keep_site_filter_list[i]!=5) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["freqs"]][i] >= 0.5) {
		    			ap <- sp_minor_allele[i]
                    } else {
		    			ap <- sp_major_allele[i]
		    		}
		    		return(ap)
		    	})
	    	unlist(per_sample_allele)        	
            })
        } else {
	        positions <- which(keep_site_filter_list==5)
	        if (!is.infinite(max_sites)) {positions <- positions[1:max_sites]}
            allele_list <- lapply(positions, function(i) {
		    	per_sample_allele <- lapply(names(site_filters), function(sample) {
		    		if (site_filters[[sample]][["flag"]][i]!=3) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["depths"]][i] == 0) {
		    			ap <- "-"
		    		} else if (site_filters[[sample]][["freqs"]][i] == -1) {
		    		    ap <- "-"	
		    		} else if (site_filters[[sample]][["freqs"]][i] >= 0.5) {
		    			ap <- sp_minor_allele[i]
		    		} else {
		    			ap <- sp_major_allele[i]
		    		}
		    		return(ap)
		    	})
	    	unlist(per_sample_allele)        	
            })
        }

        if (return_mat) {
        	mat <- do.call(cbind, allele_list)
        	row.names(mat) <- names(site_filters)
        	colnames(mat) <- row.names(SPECIES[["info"]])[positions]
	        return(mat) 
        }
        
        allele_list <- apply(do.call(cbind, allele_list), 1, function(x) paste0(x, collapse="")) |>
        setNames(names(site_filters))
        faName <- paste0(sp,"_consensus.fasta")
        if (output_seq) {
            qqcat("  Outputting consensus sequence to @{faName}\n")    	
        }
        
        fileConn<-file(faName, open="w")
        for (sample in names(site_filters)) {
			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
			cat(allele_list[sample], file=fileConn, append=TRUE, sep="\n")
		}
		close(fileConn)


		tre <- read.phyDat(faName, format = "fasta")
		stana@fastaList[[sp]] <- tre
		if (tree) {
			dm <- dist.ml(tre, "F81")
			tre <- NJ(dm)
			stana@treeList[[sp]] <- tre
			if (!is.null(cl)) {
			    tre <- groupOTU(tre, cl)
			    tp <- ggtree(tre, aes(color=.data$group),
		               layout='circular',branch.length = "none") + # Return cladogram by default
			           geom_tippoint(size=3) + ggtitle(sp)+scale_color_manual(values=stana@colors)
			    stana@treePlotList[[sp]] <- tp
			}
		}
	}
	stana
}
