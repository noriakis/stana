
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
#' @import NMF
#' @export
NMF <- function(stana, species, rank=3, target="KO", seed=53, method="snmf/r",
    deleteZeroDepth=FALSE, beta=0.01, estimate=FALSE, estimate_range=1:6, nnlm_flag=FALSE,
    nnlm_args=list()) {
	
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
    if (target=="KO") {
        mat <- stana@kos[[species]]
    } else if (target=="genes") {
        mat <- stana@genes[[species]]
    } else if (target=="snps") {
        mat <- stana@snps[[species]]
        if (deleteZeroDepth) {
           mat <- mat[rowSums(mat == -1)==0,]
           qqcat("After filtering `-1`, position numbers: @{dim(mat)[1]}\n")
        } else {
           mat[ mat == -1 ] <- NA
           nnlm_flag <- TRUE        	
        }
    }
    cat_subtle("# Original features: ", dim(mat)[1], "\n", sep="")
    cat_subtle("# Original samples: ", dim(mat)[2], "\n", sep="")
    if (!nnlm_flag) {
	    mat <- mat[apply(mat, 1, function(x) sum(x)!=0),]
	    mat <- mat[,apply(mat, 2, function(x) sum(x)!=0)]

	    cat_subtle("# Filtered features:", dim(mat)[1], "\n", sep="")
	    cat_subtle("# Filtered samples:", dim(mat)[2], "\n", sep="")
    }

    ## Test multiple ranks
    if (estimate) {
    	if (nnlm_flag) {
    		cat_subtle("# NNLM flag enabled, the cross-validation error matrix only will be returned.\n")
    		## Following the vignette procedures in NNLM
			## Comparing the MSE and randomly assigned missing values
    		A <- as.matrix(mat)
			already_na <- which(is.na(A))
			allind <- seq_len(length(A))
			newind <- allind[!(allind %in% already_na)]
			ind <- sample(newind, 0.1*length(newind));
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
	    	test <- nmfEstimateRank(as.matrix(mat),
	    		range=estimate_range, method=method)
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
			cat("Chosen rank:", rank, "\n")    		
    	}
    }
    
    cat("Rank", rank, "\n")
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

    cat("Present feature per strain:", apply(basisMat!=0, 2, function(x) sum(x)), "\n")
    stana@NMF[[species]] <- res
    return(stana)
}

#' plotStackedBarPlot
#' 
#' plot the stacked bar plot of NMF results
#' 
#' @param stana stana boject
#' @param sp species
#' @return ggplot
plotStackedBarPlot <- function(stana, sp, by="NMF") {
	if (is.null(stana@NMF[[sp]]) & is.null(stana@coefMat[[sp]])) {
		stop("NMF results or coefficient matrix should be set")
	}
	if (by=="NMF") {
    	res <- stana@NMF[[sp]]
        coefMat <- coef(res)	
	} else {
		coefMat <- stana@coefMat[[sp]]
	}
    relab <- apply(coefMat, 2, function(x) x / sum(x))
	stb <- data.frame(t(relab))
    stb[["sample"]] <- row.names(stb)
    if (length(stana@cl)!=0) {
        cols <- as.character(listToNV(stana@cl)[stb[["sample"]]])
        cc <- stana@colors %>% setNames(unique(cols))
        x_cols <- cc[cols]
    }
    melted <- reshape2::melt(stb)
    print(melted)
	ggplot(melted, aes(fill=variable, y=value, x=sample)) + 
	    geom_bar(position="fill", stat="identity")+
	    cowplot::theme_cowplot()+
	    theme(axis.text.x = element_text(angle=90, colour=x_cols))
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
        geom_boxplot()+
        facet_wrap(.~name)+
        cowplot::theme_cowplot()
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
		dat <- mat
	} else {
	  dat <- stana@NMF[[species]]
	  dat <- basis(dat)		
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
  # colnames(pathdf) <- as.character(paste0("factor",seq_len(ncol(pathdf))))
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
