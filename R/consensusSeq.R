#'
#' consensusSeq
#'
#' Output consensus sequences from merged SNV output.
#' Optionally, return phylogenetic tree inferred by `phangorn`.
#' If specified cluster of samples, additionally returns plot by `ggtree`.
#'
#' @param midas_merge_dir merge directory
#' @param species species vectors
#' @export
#'
consensusSeq <- function(target="MIDAS1", ...){
	if (target=="MIDAS2") {
		consensusSeqMIDAS2(...)
	} else if (target=="MIDAS1"){
		consensusSeqMIDAS1(...)
	} else {
		stop("please specify MIDAS1 or MIDAS2.")
	}
}

#' consensusSeqMIDAS1
#' @export
consensusSeqMIDAS1 <- function(
	midas_merge_dir,
	species,
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
    max_samples=Inf,
    keep_samples=NULL,
    exclude_samples=NULL,
    rand_samples=NULL) {
    ## site-list is currently not supported.
	files <- c("freq","depth","info","summary")
	retList <- list()
	for (sp in species) {
		qqcat("Beginning @{sp}\n")
		SPECIES <- list()
		for (file in files) {
			filePath <- paste0(midas_merge_dir,"/",sp,"/snps_",file,".txt")
			SPECIES[[file]] <- read.table(filePath, sep="\t", header=1)
		}
		siteNum <- dim(SPECIES[["freq"]])[1]
		qqcat("  Site number: @{siteNum}\n")
		SAMPLES <- list()
        for (i in seq_len(nrow(SPECIES[["summary"]]))) {
        	info <- SPECIES[["summary"]][i,]
        	if (info$fraction_covered < fract_cov) {
        		continue
        	} else if (info$mean_coverage < mean_depth) {
        		continue
        	} else {
	        	SAMPLES[[info$sample_id]] <- list(mean_depth=info$mean_coverage,
	        		fract_cov=info$fraction_covered)
        	}
        }
        for (sample in names(SAMPLES)) {
        	SAMPLES[[sample]][["freq"]] <- as.numeric(SPECIES[["freq"]][sample][,1])
        	SAMPLES[[sample]][["depth"]] <- as.numeric(SPECIES[["depth"]][sample][,1])
        	## Append sample-wise filter
	        SAMPLES[[sample]][["filter"]][["site_depth"]] <- 
	            SAMPLES[[sample]][["depth"]] >=  site_depth
	        SAMPLES[[sample]][["filter"]][["depth_ratio"]] <- 
	            (SAMPLES[[sample]][["depth"]] / SAMPLES[[sample]]$mean_depth) <= site_ratio
	        SAMPLES[[sample]][["filter"]][["allele_support"]] <- 
	            apply(cbind(SAMPLES[[sample]][["freq"]], 1-SAMPLES[[sample]][["freq"]]), 1, max) >= allele_support
	    }
        SITEFILTERS <- list()
        retainedSites <- 0
        for (i in seq_len(nrow(SPECIES[["freq"]]))) {
        	if (retainedSites >= max_sites) {break}
        	keepSamples <- c()
        	pooledMaf <- c()
        	for (sample in names(SAMPLES)) {
        		keep <- sum(SAMPLES[[sample]][["filter"]][["site_depth"]][i],
        		SAMPLES[[sample]][["filter"]][["depth_ratio"]][i],
        		SAMPLES[[sample]][["filter"]][["allele_support"]][i])
        		if (keep==3) {
        			keepSamples <- c(keepSamples, sample)
        			## Weight function is currently not supported.
        			pooledMaf <- c(pooledMaf, SAMPLES[[sample]][["freq"]][i])
        		}
        	}
			if(is.null(keepSamples)) {pooledMaf <- 0}
        	SITEFILTERS[[as.character(i)]][["prevalance"]] <- (length(keepSamples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev))
        	SITEFILTERS[[as.character(i)]][["pooled_maf"]] <- mean(pooledMaf) >= site_maf
	        SITEFILTERS[[as.character(i)]][["ref_allele"]] <- SPECIES[["info"]][i,]$ref_allele %in% c("A","T","G","C")
	        SITEFILTERS[[as.character(i)]][["site_type"]] <- SPECIES[["info"]][i,]$site_type %in% c("1D","2D","3D","4D")
	        SITEFILTERS[[as.character(i)]][["locus_type"]] <- SPECIES[["info"]][i,]$locus_type %in% c("CDS")

	        keepSite <- 
	        SITEFILTERS[[as.character(i)]][["prevalance"]]+
	        SITEFILTERS[[as.character(i)]][["pooled_maf"]]+
	        SITEFILTERS[[as.character(i)]][["ref_allele"]]+
	        SITEFILTERS[[as.character(i)]][["site_type"]]+
	        SITEFILTERS[[as.character(i)]][["locus_type"]]

	        if (keepSite==5) {
	        	retainedSites <- retainedSites + 1
	        	for (sample in names(SAMPLES)){
	        		if (!sample %in% keepSamples) {
	        			ap <- "-"
	        		} else if (SAMPLES[[sample]][["depth"]][i]==0) {
	        			ap <- "-"
	        		} else if (SAMPLES[[sample]][["freq"]][i]>=0.5) {
	        			ap <- SPECIES[["info"]][i,]$minor_allele
	        		} else {
	        			ap <- SPECIES[["info"]][i,]$major_allele
	        		}
	        		SAMPLES[[sample]][["consensus"]] <- 
	        		paste0(SAMPLES[[sample]][["consensus"]],ap)
	        	}
	        }
        }
        faName <- paste0(sp,"_consensus.fasta")
        qqcat("  Outputting consensus sequence to @{faName}\n")
        fileConn<-file(faName, open="w")
        for (sample in names(SAMPLES)) {
			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
			cat(SAMPLES[[sample]][["consensus"]], file=fileConn, append=TRUE, sep="\n")
		}
		close(fileConn)
		if (tree) {
			tre <- read.phyDat(faName, format = "fasta")
			dm <- dist.ml(tre, "F81")
			tre <- NJ(dm)
			retList[["tree"]][[sp]] <- tre
			if (!is.null(cl)) {
			    tre <- groupOTU(tre, cl)
			    tp <- ggtree(tre, aes(color=group),
		               layout='circular',branch.length = "none") + # Return cladogram by default
			           geom_tippoint(size=3) + ggtitle(sp)
			    retList[["plotTree"]][[sp]] <- tp
			}
		}
	}
	if (tree) {
		return(retList)
	}
}


#' consensusSeqMIDAS2
#' @export
consensusSeqMIDAS2 <- function(
	midas_merge_dir,
	species,
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
    max_samples=Inf,
    keep_samples=NULL,
    exclude_samples=NULL,
    rand_samples=NULL) {
    ## site-list is currently not supported.
	files <- c("freqs","depth","info")
	retList <- list()
	for (sp in species) {
		qqcat("Beginning @{sp}\n")
		SPECIES <- list()
		for (file in files) {
			filePath <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_",file,".tsv.lz4")
            filePathUn <- gsub(".lz4","",filePath)
	        system2("lz4", args=c("-d","-f",
	                                paste0(getwd(),"/",filePath),
	                                paste0(getwd(),"/",filePathUn)),
	                  stdout=FALSE, stderr=FALSE)
			SPECIES[[file]] <- read.table(filePathUn, row.names=1, header=1)
	        unlink(paste0(getwd(),"/",filePathUn))
		}

		siteNum <- dim(SPECIES[["freqs"]])[1]
		qqcat("  Site number: @{siteNum}\n")

		filePath <- paste0(midas_merge_dir,"/snps/snps_summary.tsv")
		snpsSummary <- read.table(filePath, header=1)
        SPECIES[["summary"]] <- subset(snpsSummary, snpsSummary$species_id==sp)

		SAMPLES <- list()
        for (i in seq_len(nrow(SPECIES[["summary"]]))) {
        	info <- SPECIES[["summary"]][i,]
        	if (info$fraction_covered < fract_cov) {
        		continue
        	} else if (info$mean_coverage < mean_depth) {
        		continue
        	} else {
	        	SAMPLES[[info$sample_name]] <- list(mean_depth=info$mean_coverage,
	        		fract_cov=info$fraction_covered)
        	}
        }
        for (sample in names(SAMPLES)) {
        	SAMPLES[[sample]][["freqs"]] <- as.numeric(SPECIES[["freqs"]][sample][,1])
        	SAMPLES[[sample]][["depth"]] <- as.numeric(SPECIES[["depth"]][sample][,1])
        	## Append sample-wise filter
	        SAMPLES[[sample]][["filter"]][["site_depth"]] <- 
	            SAMPLES[[sample]][["depth"]] >=  site_depth
	        SAMPLES[[sample]][["filter"]][["depth_ratio"]] <- 
	            (SAMPLES[[sample]][["depth"]] / SAMPLES[[sample]]$mean_depth) <= site_ratio
	        SAMPLES[[sample]][["filter"]][["allele_support"]] <- 
	            apply(cbind(SAMPLES[[sample]][["freqs"]], 1-SAMPLES[[sample]][["freqs"]]), 1, max) >= allele_support
	    }
        SITEFILTERS <- list()
        retainedSites <- 0
        for (i in seq_len(nrow(SPECIES[["freqs"]]))) {
        	if (retainedSites >= max_sites) {break}
        	keepSamples <- c()
        	pooledMaf <- c()
        	for (sample in names(SAMPLES)) {
        		keep <- sum(SAMPLES[[sample]][["filter"]][["site_depth"]][i],
        		SAMPLES[[sample]][["filter"]][["depth_ratio"]][i],
        		SAMPLES[[sample]][["filter"]][["allele_support"]][i])
        		if (keep==3) {
        			keepSamples <- c(keepSamples, sample)
        			## Weight function is currently not supported.
        			pooledMaf <- c(pooledMaf, SAMPLES[[sample]][["freqs"]][i])
        		}
        	}
			if(is.null(keepSamples)) {pooledMaf <- 0}
        	SITEFILTERS[[as.character(i)]][["prevalance"]] <- (length(keepSamples) / length(names(SAMPLES))) >= max(c(1e-6, site_prev))
        	SITEFILTERS[[as.character(i)]][["pooled_maf"]] <- mean(pooledMaf) >= site_maf
	        # SITEFILTERS[[as.character(i)]][["ref_allele"]] <- SPECIES[["info"]][i,]$ref_allele %in% c("A","T","G","C")
	        SITEFILTERS[[as.character(i)]][["site_type"]] <- SPECIES[["info"]][i,]$site_type %in% c("1D","2D","3D","4D")
	        SITEFILTERS[[as.character(i)]][["locus_type"]] <- SPECIES[["info"]][i,]$locus_type %in% c("CDS")

	        keepSite <- 
	        SITEFILTERS[[as.character(i)]][["prevalance"]]+
	        SITEFILTERS[[as.character(i)]][["pooled_maf"]]+
	        SITEFILTERS[[as.character(i)]][["site_type"]]+
	        SITEFILTERS[[as.character(i)]][["locus_type"]]

	        if (keepSite==4) {
	        	retainedSites <- retainedSites + 1
	        	for (sample in names(SAMPLES)){
	        		if (!sample %in% keepSamples) {
	        			ap <- "-"
	        		} else if (SAMPLES[[sample]][["depth"]][i]==0) {
	        			ap <- "-"
	        		} else if (SAMPLES[[sample]][["freqs"]][i]>=0.5) {
	        			ap <- SPECIES[["info"]][i,]$minor_allele
	        		} else {
	        			ap <- SPECIES[["info"]][i,]$major_allele
	        		}
	        		SAMPLES[[sample]][["consensus"]] <- 
	        		paste0(SAMPLES[[sample]][["consensus"]],ap)
	        	}
	        }
        }
        faName <- paste0(sp,"_consensus.fasta")
        qqcat("  Outputting consensus sequence to @{faName}\n")
        fileConn<-file(faName, open="w")
        for (sample in names(SAMPLES)) {
			cat(paste0(">",sample), file=fileConn, append=TRUE, sep="\n")
			cat(SAMPLES[[sample]][["consensus"]], file=fileConn, append=TRUE, sep="\n")
		}
		close(fileConn)
		if (tree) {
			tre <- read.phyDat(faName, format = "fasta")
			dm <- dist.ml(tre, "F81")
			tre <- NJ(dm)
			retList[["tree"]][[sp]] <- tre
			if (!is.null(cl)) {
			    tre <- groupOTU(tre, cl)
			    tp <- ggtree(tre, aes(color=group),
		               layout='circular',branch.length = "none") + # Return cladogram by default
			           geom_tippoint(size=3) + ggtitle(sp)
			    retList[["plotTree"]][[sp]] <- tp
			}
		}
	}
	if (tree) {
		return(retList)
	}
}