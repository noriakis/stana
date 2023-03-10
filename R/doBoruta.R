#' doBoruta
#' 
#' feature selection of classification between category
#' using snv frequency and gene matrix
#' 
#' @param stana stana object
#' @param sp candidate species ID
#' @param cl named list of category
#' @param target default to snps
#' @param mat if target is not snps, provide preprocessed gene matrix
#' otherwise the raw gene matrix is used.
#' @param whichToCount which to show on box plots,
#' pass to checkPATRIC() (in MIDAS1)
#' @param argList passed to Boruta()
#' @return Boruta object
#' @import ggplot2
#' @export
doBoruta <- function(stana, sp, cl,
  target="snps", mat=NULL, whichToCount="ec_description", argList=list()) {
  ret <- list()
  if (target!="snps" & is.null(mat)){
    qqcat("If needed, please provide preprocessed matrix of genes\n")
    filtDf <- stana@genes[[sp]]
  }
  if (target=="snps"){
    filtDf <- stana@snps[[sp]]
  } else {
    qqcat("Proceeding with provided matrix\n")
    filtDf <- mat
  }
  transDf <- data.frame(t(filtDf),
                        check.names=F)
  gr <- NULL
  for (cn in rownames(transDf)){
    for (clm in seq_along(cl)){
      if (cn %in% cl[[clm]]) {
        gr <- c(gr, names(cl)[clm])
      }
    }
  }
  
  transDf$group <- as.factor(gr)
  qqcat("Performing Boruta\n")
  argList[["formula"]] <- group ~ .
  argList[["data"]] <- transDf
  rf <- do.call("Boruta", argList)
  bor <- Boruta::TentativeRoughFix(rf)
  ret[["boruta"]] <- bor
  
  dec <- bor$finalDecision
  confirmedIDs <- names(dec[dec=="Confirmed"])
  if (length(confirmedIDs)==0) {stop("No feature selected")}
  if (target=="snps" & stana@type=="MIDAS1") {
    info <- read.table(paste0(stana@mergeDir,"/",cand,"/snps_info.txt"),
                       sep="\t", header=1)
    whichGene <- info[ info$site_id %in% confirmedIDs, ]$gene_id
    whichGene <- whichGene[ !is.na(whichGene) ]
  } else if (target!="snps"){
    whichGene <- confirmedIDs
  } else {
    qqcat("currently not supported\n")
    return(ret)
  }
  
  if (target=="copynum") {
    pl <- transDf[,c(whichGene,"group")]
    pll <- tidyr::pivot_longer(pl, seq_len(length(whichGene)))
    
    if (stana@type=="MIDAS1"){
      pat <- checkPATRIC(list(important=whichGene),
                         whichToCount=whichToCount)
      patnm <- unlist(lapply(strsplit(pat$important$DF$patric_id,"\\|"),"[",2))
      desc <- list()
      for (i in colnames(pl)) {
        tmp <- pat$important$DF[patnm %in% i,] 
        if (dim(tmp)[1]!=0){
          desc[[i]] <-
            paste0(unique(tmp$ec_description),
                   collapse=" / ")
        } else {
          desc[[i]] <- NULL
        }
      }
      newname <- NULL
      for (i in pll$name) {
        if ( !is.null(desc[[i]]) ) {
          newname <- c(newname, stringr::str_wrap(
            paste0(i, "\n(", desc[[i]],")"), 20))
        } else {
          newname <- c(newname, NA)
        }
      }
      pll$newname <- newname
      pllna <- pll[!is.na(pll$newname),]
      box<-ggplot(pllna, aes(x=.data$group, y=.data$value,
                             fill=.data$group, color=.data$group))+
        geom_jitter()+geom_boxplot(alpha=0.2)+
        facet_wrap(.~newname, nrow=1)+
        theme_minimal()
      ret[["plot"]] <- box
    } else {
      box<-ggplot(pll, aes(x=.data$group, y=.data$value,
                           fill=.data$group, color=.data$group))+
        geom_jitter()+geom_boxplot(alpha=0.2)+
        facet_wrap(.~name)+
        theme_minimal()
      ret[["plot"]] <- box
    }
  }
  return(ret)
}
