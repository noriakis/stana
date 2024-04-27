#'
#' doAdonis
#' 
#' Perform PERMANOVA on distance matrix based 
#' on various distances based on such as SNV frequency or 
#' gene matrix using `adonis2`.
#' Named list of samples are to be provided.
#' Note that this function currently performs 
#' comparison using one category.
#' for complex design, one should do it manually.
#' 
#' @param stana stana object
#' @param specs species to be examined
#' @param cl named list of samples
#' @param target tree, snps, genes, fasta
#' @param formula pass to adonis2, specify distance matrix as d.
#' @param distMethod distance method passed to dist() (default, manhattan)
#' @param maj major allele distance
#' @param AAfunc if choose `fasta`, provide function for calculating distance
#' @param AAargs provided to `AAfunc`
#' @param argList parameters passed to adonis2
#' @param deleteZeroDepth delete zero depth snvs (denoted as `-1`)
#' Otherwise the cell is treated as NA
#' @param distArg passed to `dist()`
#' @param useMeta use metadata slot, not cl slot
#' @param pcoa if TRUE, performs ape::pcoa on dist and plot them
#' @importFrom vegan adonis2
#' @importFrom stats as.formula dist
#' @importFrom utils read.table
#' @export
doAdonis <- function(stana, specs, cl=NULL,
    target="snps", formula=NULL,
    distMethod="manhattan", pcoa=FALSE,
    AAfunc=dist.ml, AAargs=list(),
    maj=FALSE, deleteZeroDepth=FALSE,
    useMeta=FALSE,
    argList=list(), distArg=list()) {
      if (is.null(cl)) {cl <- stana@cl}
      if (pcoa) {pcoaList<-list()}
      for (sp in specs){
        cat_subtle("# Performing adonis in ", sp, " target is ", target, "\n", sep="")
        if (target=="snps") {
            snps <- stana@snps[[sp]]
			if (!is.null(stana@includeSNVID[[sp]])) {
				cat_subtle("# The set SNV ID information (", length(stana@includeSNVID[[sp]]), ") is used.\n")
				snps <- snps[stana@includeSNVID[[sp]], ]
			}
            if (deleteZeroDepth) {
              snps <- snps[rowSums(snps==-1)==0,]
              cat_subtle("After filtering `-1`, position numbers: ", dim(snps)[1], "\n", sep="")
            }
            snps[snps == -1] <- NA
        } else if (target=="tree") {
            if (!is.null(stana@treeList[[sp]])) {
              tre <- stana@treeList[[sp]]
              ## cophenetic distance
              d <- as.dist(ape::cophenetic.phylo(tre))
              sn <- attr(d, "Labels")
            } else {
              stop("No tree found in stana@treeList")
            }
        } else if (target=="fasta") {
          AAargs[["x"]] <- stana@fastaList[[sp]]
          d <- do.call(AAfunc, AAargs)
          sn <- attr(d, "Labels")
        } else {
            snps <- stana@genes[[sp]]          
        }

        if (!(target %in% c("tree","fasta"))) {
          if (maj & target=="snps" & stana@type=="MIDAS1") {
              chk <- read.table(paste0(stana@mergeDir,
                "/",sp,"/snps_info.txt"),
              header=1)
              qqcat("  before filtering: @{dim(chk)[1]}\n")
              calcDif <- function(x){
                if (x[14]=="bi") {
                    maj <- paste0("count_",tolower(x[5]))
                    min <- paste0("count_",tolower(x[6]))
                    fre <- as.numeric(x[8:11]) / sum(as.numeric(x[8:11]))
                    names(fre) <- names(x[8:11])
                    fre[maj]-fre[min] > 0.6
                } else {
                  FALSE
                }
              }
            filtIDs <- chk[apply(chk, 1, function(x) calcDif(x)),]$site_id
            qqcat("  after filtering: @{length(filtIDs)}\n")

            inc <- intersect(filtIDs, row.names(snps))
            snps <- snps[inc,]
          }
          snps <- snps[,intersect(colnames(snps),as.character(unlist(cl)))]
          distArg[["x"]] <- t(snps)
          distArg[["method"]] <- distMethod
          d <- do.call(dist, distArg)
          sn <- attr(d, "Labels")
          # gr <- as.character(listToNV(stana@cl)[sn])
        	met <- data.frame(row.names=sn)
        	met[["group"]] <- listToNV(stana@cl)[row.names(met)]
          # gr <- NULL
          # for (cn in sn){
          # 	if (cn %in% unlist(cl)) {
	      #       for (clm in seq_along(cl)){
	      #         if (cn %in% cl[[clm]]) {
	      #           gr <- c(gr, names(cl)[clm])
	      #         }
	      #       }          		
          # 	} else {
          # 		gr <- c(gr, NA)
          # 	}
          # }       
        } else {
        	met <- data.frame(row.names=sn)
        	met[["group"]] <- listToNV(stana@cl)[row.names(met)]
          # gr <- NULL
          # for (cn in sn){
          # 	if (cn %in% unlist(cl)) {
	      #       for (clm in seq_along(cl)){
	      #         if (cn %in% cl[[clm]]) {
	      #           gr <- c(gr, names(cl)[clm])
	      #         }
	      #       }          		
          # 	} else {
          # 		gr <- c(gr, NA)
          # 	}
          # }
          # d <- as.dist(d)         
        }
        if (is.null(formula)){
            formulaPass <- as.formula("d ~ .")
            pr <- TRUE
        } else {
            formulaPass <- as.formula(formula)
            pr <- FALSE
        }
        argList[["formula"]] <- formulaPass
        if (useMeta) {
          ## Assuming all the samples are in the metadata
          argList[["data"]] <- stana@meta[sn,]
          if (is.null(formula)) {
            stop("Please specify formula when using metadata slot")
          }
        } else {
          argList[["data"]] <- met
        }

        adores <- do.call("adonis2", argList)
        
        if (pr) {
            pr <- adores$`Pr(>F)`
            pr <- pr[!is.na(pr)]
            r2 <- adores$R2[1]
            cat_subtle("#  F: ",adores$F[1],", R2: ",r2,", Pr: ",pr,"\n",sep="")
        } else {
          cat_subtle("# Printing raw adonis results ...\n")
          print(adores)
        }
        stana@adonisList[[sp]] <- adores
        
        if (pcoa) {
        	pco <- ape::pcoa(d)
        	dat <- data.frame(pco$vector[, 1:2])
        	relei <- pco$values$Relative_eig
        	dat[["group"]] <- listToNV(stana@cl)[row.names(dat)]
        	pcoaList[[sp]] <- ggplot(dat, aes(x=dat[,1], y=dat[,2])) +
            	geom_point(aes(color=group))+
            	xlab(paste0("Axis 1 (", round(relei[1], 3),")"))+
            	ylab(paste0("Axis 2 (", round(relei[2], 3),")"))+
            	scale_color_manual(values=stana@colors)+
            	cowplot::theme_cowplot()
        }
      }
      if (pcoa) {
          plot(patchwork::wrap_plots(pcoaList))  	
      }
      stana
}