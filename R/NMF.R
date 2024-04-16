
#' NMF
#'
#' decompose SNV, gene content, or gene family abundance to
#' factor x sample and factor x feature matrix.
#' 
#' @param stana stana object
#' @param species candidate species ID
#' @param rank rank of NMF
#' @param target KO, gene or snv, default to KO
#' @param seed random seed
#' @param method NMF method, default to snmf/r
#' @param beta argument to be passed to NMF function
#' @param deleteZeroDepth when snv matrix is used, this option filters the
#' positions with zero-depth (indicated by -1)
#' @param estimate estimate rank
#' @param estimate_range range of ranks for the estimation
#' @param nnlm_flag if TRUE, use NNLM package which can handle missing values in the data.
#' The different approach for estimating k is taken when estimate is TRUE.
#' When estimate is TRUE, only the error matrix is returned.
#' @param nnlm_args arguments passed to NNLM functions
#' @param nnlm_na_perc NA frequency when the cross-validation is performed.
#' @param tss perform total sum scaling to the matrix
#' @param remove_na_perc remove features with too many NA
#' (features with NA below {remove_na_perc} * overall sample numbers will be retained)
#' @importFrom NMF nmf nmfEstimateRank
#' @export
NMF <- function(stana, species, rank=3, target="kos", seed=53, method="snmf/r",
    deleteZeroDepth=FALSE, beta=0.01, estimate=FALSE, estimate_range=1:6, nnlm_flag=FALSE,
    nnlm_args=list(), nnlm_na_perc=0.3, tss=FALSE, remove_na_perc=NULL) {
	
	## References for the missing value handling in NMF:
	## https://github.com/scikit-learn/scikit-learn/pull/8474
	## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6945623/
	
	cat_subtle("# NMF started ", species, ", target: ", target, ", method: ", ifelse(nnlm_flag, "NNLM::nnmf", method), "\n", sep="")
    if (length(species)>1) {stop("NMF accepts only one species per run")}
    if (nnlm_flag) {
    	if (!requireNamespace("NNLM")) {
    		stop("Please install NNLM.")
    	}
    }
    if (target=="kos") {
        mat <- stana@kos[[species]]
    } else if (target=="genes") {
        mat <- stana@genes[[species]]
    } else if (target=="snps") {
    	## For SNV, zero-inflated models should be applied (with NA handling capabilities)
        mat <- stana@snps[[species]]
		if (!is.null(stana@includeSNVID[[species]])) {
			cat_subtle("# The set SNV ID information (", length(stana@includeSNVID[[species]]), ") is used.\n")
			mat <- mat[stana@includeSNVID[[species]], ]
		}
        if (deleteZeroDepth) {
           mat <- mat[rowSums(mat == -1)==0,]
           cat_subtle("# After filtering `-1`, position numbers: ",dim(mat)[1],"\n")
        } else {
           mat[ mat == -1 ] <- NA
           # if (!nnlm_flag) cat_subtle("# Changing to NNLM\n")
           # nnlm_flag <- TRUE
        }
    }
    
    nac <- as.numeric(table(is.na(mat))["TRUE"]) / (dim(mat)[1]*dim(mat)[2])
    zeroc <- as.numeric(table(mat==0)["TRUE"]) / (dim(mat)[1]*dim(mat)[2])
    
    cat_subtle("# Original features: ", dim(mat)[1], "\n", sep="")
    cat_subtle("# Original samples: ", dim(mat)[2], "\n", sep="")
    cat_subtle("# Original matrix NA: ", round(nac, 3), "\n", sep="")
    cat_subtle("# Original matrix zero: ", round(zeroc, 3), "\n", sep="")

    if (is.null(nnlm_args[["loss"]])){
    	cat_subtle("# Selecting KL loss\n")
       	nnlm_args[["loss"]] <- "mkl"
    }
    
    if (!is.null(remove_na_perc)){
    	cat_subtle("# Removing too many NA features: ", remove_na_perc * dim(mat)[2], "\n", sep="")
    	mat <- mat[apply(is.na(mat), 1, sum) < remove_na_perc * dim(mat)[2],]
    }
    
    if (tss) {
    	cat_subtle("# Performing TSS\n")
    	mat <- apply(mat, 2, function(x) x / sum(x))
    }
    if (!nnlm_flag) {
	    mat <- mat[apply(mat, 1, function(x) sum(x)!=0),]
	    mat <- mat[,apply(mat, 2, function(x) sum(x)!=0)]
    }
    cat_subtle("# Filtered features: ", dim(mat)[1], "\n", sep="")
    cat_subtle("# Filtered samples: ", dim(mat)[2], "\n", sep="")

    ## Test multiple ranks
    if (estimate) {
    	if (nnlm_flag) {
    		cat_subtle("# NNLM flag enabled, the error matrix only will be returned.\n")
    		## Following the vignette procedures in NNLM
			## Comparing the MSE and randomly assigned missing values
    		A <- as.matrix(mat)
			already_na <- which(is.na(A))
			allind <- seq_len(length(A))
			newind <- allind[!(allind %in% already_na)]
			ind <- sample(newind, nnlm_na_perc * length(newind));
			A2 <- A;
			A2[ind] <- NA;

			err <- sapply(X = estimate_range,
			              FUN = function(k) {
			                z <- nnmf(A2, k);
			                c(mean((with(z, W %*% H)[ind] - A[ind])^2), tail(z$mse, 1));
			              }
			);
			return(err)
    	} else {
    		## Following the cophenetic correlation coefficient drop procedure
            if (any(is.na(mat))) {
                ## If NA included, use ls-nmf with weight
                w <- matrix(1, nrow(mat), ncol(mat))
                w[ is.na(mat) ] <- 0
                mat[is.na(mat)] <- 123456789
                test <- nmfEstimateRank(as.matrix(mat),
                    range=estimate_range, method="ls-nmf",
                    weight=w)
            } else {
                test <- nmfEstimateRank(as.matrix(mat),
                    range=estimate_range, method=method)
            }
	    	val <- test$measures[, "cophenetic"]
	        b <- -1
			for (i in seq_along(val)) {
			  if (is.na(val[i])) {
			    next
			  } else {
			    if (val[i] > b) {
			      b <- val[i]
			    } else {
			      break
			    }
			  }
			}
		    rank <- estimate_range[i]
			cat_subtle("# Chosen rank:", rank, "\n")    		
    	}
    }
    
    cat_subtle("# Rank ", rank, "\n", sep="")
    if (nnlm_flag){
    	nnlm_args[["A"]] <- as.matrix(mat)
    	nnlm_args[["k"]] <- rank
    	res <- do.call(NNLM::nnmf, nnlm_args)
    	# res <- NNLM::nnmf(mat, rank)
    	coefMat <- res$H
    	basisMat <- res$W
    } else {
	    if (method %in% c("snmf/l", "snmf/r")) {
	        res <- NMF::nmf(mat, rank = rank, seed = seed, method=method, beta=beta)
	    } else {
	        res <- NMF::nmf(mat, rank = rank, seed = seed, method=method)
	    }
	    coefMat <- coef(res)
        basisMat <- basis(res) 	
    }

    stana@coefMat[[species]] <- data.frame(coefMat)
    ## Plot by default
    relab <- apply(coefMat, 2, function(x) x / sum(x))
    cat("Mean relative abundances:", apply(relab, 1, mean), "\n")

    cat("Present feature per factor:", apply(basisMat!=0, 2, function(x) sum(x)), "\n")
    stana@NMF[[species]] <- res
    if (estimate & !nnlm_flag) {
        return(list(stana, test))
    } else {
        return(stana)    
    }
}

#' plotStackedBarPlot
#' 
#' plot the stacked bar plot of NMF results
#' 
#' @param stana stana boject
#' @param sp species
#' @param by "NMF" or "coef"
#' @return ggplot
#' @export
plotStackedBarPlot <- function(stana, sp, by="NMF") {
	if (is.null(stana@NMF[[sp]]) & is.null(stana@coefMat[[sp]])) {
		stop("NMF results or coefficient matrix should be set")
	}
	if (by=="NMF") {
    	res <- stana@NMF[[sp]]
        coefMat <- coef(res)	
	} else if (by=="coef") {
		coefMat <- stana@coefMat[[sp]]
	} else {
		stop("NMF results or coefficient matrix should be set")
	}
    relab <- apply(coefMat, 2, function(x) x / sum(x))
	stb <- data.frame(t(relab))
	colnames(stb) <- as.character(seq_len(ncol(stb)))
    stb[["sample"]] <- row.names(stb)
    ## This does not use the stana@colors slot
    if (length(stana@cl)!=0) {
        stb[["group"]] <- as.character(listToNV(stana@cl)[stb[["sample"]]])
	    melted <- reshape2::melt(stb)
	    ## Not showing sample label, instead facet by group
		ggplot(melted, aes(fill=variable, y=value, x=sample)) + 
		    geom_col(position="fill")+
		    facet_grid(. ~ group, scale="free")+
            scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
            cowplot::theme_cowplot()+ cowplot::panel_border()+
		    theme(axis.text.x = element_blank())
    } else {
	    melted <- reshape2::melt(stb)
		ggplot(melted, aes(fill=variable, y=value, x=sample)) + 
		    geom_bar(position="fill", stat="identity")+
            scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
		    cowplot::theme_cowplot()+cowplot::panel_border()+
		    theme(axis.text.x = element_text(angle=90))    	
    }

}

#' alphaDiversityWithinSpecies
#' @param stana stana object
#' @param species species
#' @param method method for vegan::diversity
#' if `spc`, factor count will be returned.
#' @param rank if NMF is not performed, this performs the NMF beforehand.
#' rank can be specified here.
#' @export
alphaDiversityWithinSpecies <- function(stana, species, method="shannon", rank=5) {
    if (is.null(stana@NMF[[species]])) {
        stana <- NMF(stana, species, rank=rank)
    }
    res <- stana@NMF[[species]]
    H <- coef(res)
    if (method=="spc") {
        div <- apply(coef(stana@NMF[[1]])==0, 2, function(x) sum(x))
    } else {
        div <- vegan::diversity(t(H), index=method)
    }
    
    if (!is.null(stana@cl)) {
        nm <- listToNV(stana@cl)
        div <- data.frame(div)
        colnames(div) <- "alpha_diversity"
        div[["group"]] <- nm[row.names(div)]
    }
    return(div)
}

#' plotAbundanceWithinSpecies
#' 
#' plot abundances using factor to sample matrix produced by NMF.
#' 
#' @param stana stana object
#' @param species species ID
#' @param tss perform total sum scaling
#' @param return_data return only the data, not plot
#' @param by NMF or coef matrix set to `coefMat` slot
#' @export
plotAbundanceWithinSpecies <- function(stana, species, tss=TRUE, return_data=FALSE, by="NMF") {
	if (by=="NMF") {
	    if (is.null(stana@NMF[[species]])) {
	        stana <- NMF(stana, species)
    	}
    	res <- stana@NMF[[species]]
        H <- coef(res)	
	} else if (by=="coef") {
		H <- stana@coefMat[[species]]
	} else {
		stop("NMF or coef should be specified in `by`")
	}

    if (tss) {
        H <- apply(H, 2, function(x) x / sum(x))
    }
    H <- data.frame(t(H))
    if (!is.null(stana@cl)) {
        nm <- listToNV(stana@cl)
        H[["group"]] <- nm[row.names(H)]
    }
    colnames(H) <- c(as.character(seq_len(ncol(H)-1)),"group")

    if (return_data) {
        return(H)
    }
    H %>% tidyr::pivot_longer(1:(ncol(H)-1)) %>%
        ggplot(aes(x=group, y=value))+
        geom_boxplot(aes(fill=group), alpha=0.5)+
        facet_wrap(.~name)+
        cowplot::theme_cowplot()+
        cowplot::panel_border()+
        scale_fill_manual(values=stana@colors)
}


#' pathwayWithFactor
#' convert KO matrix per factor obtained by NMF function
#' to pathway to factor matrix by summing the KOs in the pathway.
#' @param stana stana boject
#' @param species species ID
#' @param tss perform total sum scaling to the resulting matrix
#' @param change_name change pathway names to description
#' @param summarize summarizing function, default to base::sum
#' @param mat other matrix than the basis of NMF
#' @export
pathwayWithFactor <- function(stana, species, tss=FALSE, change_name=FALSE,
	summarize=sum, mat=NULL) {
	if (!is.null(mat)) {
		use_name <- TRUE
		dat <- mat
	} else {
	  dat <- stana@NMF[[species]]
	  dat <- basis(dat)		
	  use_name <- FALSE
	}

  bfc <- BiocFileCache()
  url <- bfcrpath(bfc,"https://rest.kegg.jp/link/ko/pathway")

  summed <- data.frame(data.table::fread(url, header=FALSE))
  summed <- summed[grepl("ko", summed$V1),]
  ## No global map
  summed <- summed[!grepl("ko011", summed$V1),]
  summed <- summed[!grepl("ko012", summed$V1),]

  allpath <- unique(summed$V1)
  pathdf <- do.call(rbind, lapply(allpath, function(i) {
    tmp <- summed[summed$V1==i, ]
    int <- length(intersect(row.names(dat), tmp$V2))
    if (int>1) {
      tmpsum <- apply(dat[intersect(row.names(dat), tmp$V2),], 2, summarize)
      return(c(i, tmpsum))
    } else if (int==1) {
      return(c(i, dat[intersect(row.names(dat), tmp$V2),]))
    } else {
      return(NULL)
    }
  })) %>% data.frame()
  row.names(pathdf) <- pathdf[,1]
  pathdf[,1] <- NULL
  if (!use_name) {
     colnames(pathdf) <- as.character(paste0("factor",seq_len(ncol(pathdf))))	
  }
  pathdf <- dplyr::mutate_all(pathdf, as.numeric)  
  if (tss) {
    pathdf <- apply(pathdf, 2, function(x) x / sum(x))
  }
  if (change_name) {
  	  url2 <- bfcrpath(bfc,"https://rest.kegg.jp/list/pathway")
      namec <- data.frame(data.table::fread(url2, header=FALSE))
      namec <- namec$V2 %>% setNames(namec$V1)
      row.names(pathdf) <- namec[gsub("path:ko","map",row.names(pathdf))]
  }
  return(pathdf)
}
