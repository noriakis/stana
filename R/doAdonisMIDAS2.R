#'
#' doAdonisMIDAS2
#' 
#' Perform PERMANOVA on distance matrix based on SNV frequency.
#' Named list of samples are to be provided.
#' Note that this function performs comparison using one category.
#' for complex design, one should do it manually.
#' The location of `lz4` binary must be in PATH.
#' 
#' @param midas_merge_dir output directory of midas2 merge
#' @param specs species to be examined
#' @param cl named list of samples
#' @param taxtbl data.frame of 6-digit MIDAS species IDs 
#' as row.names, and have `GTDB species` column
#' @param formula pass to adonis2, specify distance matrix as d.
#' @param distMethods distance method passed to dist() (default, manhattan)
#' @param ... parameters passed to adonis2
#' @import vegan
#' @export
doAdonisMIDAS2 <- function(midas_merge_dir,
                           specs, cl, target="snps",
                           taxtbl=NULL, formula=NULL,
                           distMethod="manhattan", ...) {
  res <- list()
  for (sp in specs){
    
    qqcat("Performing adonis in @{sp}\n")
    if (!is.null(taxtbl)){
      spnm <- taxtbl[sp,]$`GTDB species`
      qqcat("  @{spnm}\n")
    }
    if (target=="snps") {
        cnc <- paste0(midas_merge_dir,"/snps/",sp,"/",sp,".snps_freqs.tsv.lz4")
    } else if (target=="presabs") {
        cnc <- paste0(midas_merge_dir,"/genes/",sp,"/",sp,".genes_copynum.tsv.lz4")
    } else if (target=="copynum") {
        cnc <- paste0(midas_merge_dir,"/genes/",sp,"/",sp,".genes_presabs.tsv.lz4")
    } else {
        stop('"Please specify one of "snps", "presabs", or "copynum"')
    }
    cnd <- gsub(".lz4","",cnc)
    system2("lz4", args=c("-d","-f",
                          paste0(getwd(),"/",cnc),
                          paste0(getwd(),"/",cnd)),
            stdout=FALSE, stderr=FALSE)
    snps <- read.table(cnd, row.names=1, header=1)
    unlink(paste0(getwd(),"/",cnd))
    d <- dist(t(snps), method=distMethod)
    gr <- c()
    for (cn in colnames(snps)){
      for (clm in seq_along(cl)){
        if (cn %in% cl[[clm]]) {
          gr <- c(gr, names(cl)[clm])
        }
      }
    }
    if (is.null(formula)){
        formulaPass <- as.formula("d ~ gr")
        pr <- TRUE
    } else {
        formulaPass <- as.formula(formula)
        pr <- FALSE
    }
    adores <- adonis2(formulaPass, ...)
    if (pr) {
        pr <- adores$`Pr(>F)`
        pr <- pr[!is.na(pr)]
        r2 <- adores$R2[1]
        qqcat("  R2: @{r2}, Pr: @{pr}\n")
    }
    res[[sp]][["adonis"]] <- adores
  }
  res
}
