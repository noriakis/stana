
#' @noRd
alleleStat_metaSNV <- function(stana, sp, cl, deleteZeroDepth) {
	statList <- list()
	freqDf <- stana@snps[[sp]]
	qqcat("Overall, @{dim(freqDf)[1]} SNVs\n")
	if (deleteZeroDepth) {
		freqDf <- freqDf[rowSums(freqDf==-1)==0,]
		qqcat("  After removal of -1, @{dim(freqDf)[1]} SNVs\n")
	}
	trans_raw <- row.names(freqDf) |> strsplit(":") |> vapply("[", 4, FUN.VALUE="character")
	statList[["transTable"]] <- table(trans_raw)
	if (length(cl)!=0) {
		freqs <- NULL
		wholeFreqs <- NULL
		for (cln in names(cl)) {
		  tmp <- freqDf[,intersect(colnames(freqDf),cl[[cln]])]
		  meanFreq <- apply(tmp, 1, mean)
		  names(meanFreq) <- trans_raw
		  meanTrans <- data.frame(tapply(meanFreq, names(meanFreq), mean))
		  ord <- row.names(meanTrans)
		  freqs <- cbind(freqs, meanTrans[ord,])
		}
		freqs <- data.frame(freqs) |>
		  `colnames<-`(names(cl)) |>
		  `row.names<-`(ord)
	    statList[["meanMafPerGroup"]] <- freqs		
	}
	return(statList)
}



#' alleleStat
#' 
#' Output statistics of MAF, using site-by-sample 
#' minor allele frequency table output by MIDAS1 and MIDAS2.
#' 
#' @param stana stana object
#' @param sp candidate species
#' @param cl named list of clusters of samples
#' @param base "maj","ref", or "majref"
#' @param deleteZeroDepth delete snvs with zero depth
#' @export
#' 

alleleStat <- function(stana, sp, cl=NULL, base="maj",
	deleteZeroDepth=FALSE) {
	if (length(sp)!=1) {stop("please provide one species")}
	if (is.null(cl)) {cl <- stana@cl}
	if (stana@type=="metaSNV") {return(alleleStat_metaSNV(stana, sp, cl, deleteZeroDepth))}
	statList <- list()
	merge_dir <- stana@mergeDir
	freqDf <- stana@snps[[sp]]
	if (is.null(stana@snpsInfo[[sp]])) {
		if (stana@type=="MIDAS1") {
			info <- read.table(paste0(merge_dir,"/",sp,
	                          "/snps_info.txt"),
	                   row.names=1, header=1)
		} else if (stana@type=="MIDAS2") {
			cnc <- paste0(merge_dir,"/snps/",sp,"/",sp,".snps_info.tsv.lz4")
			cnd <- gsub(".lz4","",cnc)
			system2("lz4", args=c("-d","-f",
			                      paste0(getwd(),"/",cnc),
			                      paste0(getwd(),"/",cnd)),
			        stdout=FALSE, stderr=FALSE)
			info <- read.table(cnd, row.names=1, header=1)
			unlink(paste0(getwd(),"/",cnd))
			## Append ref_allele info if available (default database)
			info$ref_allele <- sapply(strsplit(row.names(info), "\\|"), "[", 5)
		} else {
			stop("Currently, MIDAS1 and MIDAS2 is supported")
		}		
	} else {
		info <- stana@snpsInfo[[sp]]
	}
	qqcat("Overall, @{dim(info)[1]} SNVs\n")
	if (deleteZeroDepth) {
		freqDf <- freqDf[rowSums(freqDf==-1)==0,]
		qqcat("  After removal of -1, @{dim(freqDf)[1]} SNVs\n")
	}
	if (base=="ref") {
	  info$trans <- paste0(info$ref_allele,">",info$minor_allele)
	} else if (base=="maj") {
	  info$trans <- paste0(info$major_allele,">",info$minor_allele)
	} else {
	  info$trans <- paste0(info$ref_allele,":",info$major_allele,">",info$minor_allele)
	}
	statList[["transTable"]] <- table(info$trans)
	if (length(cl)!=0) {
		freqs <- NULL
		wholeFreqs <- NULL
		for (cln in names(cl)) {
		  tmp <- freqDf[,intersect(colnames(freqDf),cl[[cln]])]
		  meanFreq <- apply(tmp, 1, mean)
		  names(meanFreq) <- info[names(meanFreq),"trans"]
		  meanTrans <- data.frame(tapply(meanFreq, names(meanFreq), mean))
		  ord <- row.names(meanTrans)
		  freqs <- cbind(freqs, meanTrans[ord,])
		}
		freqs <- data.frame(freqs) |>
		  `colnames<-`(names(cl)) |>
		  `row.names<-`(ord)
	    statList[["meanMafPerGroup"]] <- freqs		
	}
	return(statList)
}