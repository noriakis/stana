#'
#' doAdonisMIDAS1
#' 
#' Perform PERMANOVA on distance matrix based on SNV frequency.
#' Named list of samples are to be provided.
#' Note that this function performs comparison using one category.
#' for complex design, one should do it manually.
#' 
#' @param midas_merge_dir output directory of midas2 merge
#' @param specs species to be examined
#' @param cl named list of samples
#' @param target snps, presabs, or copynum
#' @param formula pass to adonis2, specify distance matrix as d.
#' @param distMethod distance method passed to dist() (default, manhattan)
#' @param ... parameters passed to adonis2
#' @importFrom vegan adonis2
#' @importFrom stats as.formula dist
#' @importFrom utils read.table
#' @export
doAdonisMIDAS1 <- function(midas_merge_dir, specs, cl,
    target="snps", formula=NULL, distMethod="manhattan", ...) {
      res <- list()
      for (sp in specs){
        qqcat("Performing adonis in @{sp}\n")
        if (target=="snps") {
            snps <- read.table(paste0(midas_merge_dir,"/",sp,"/snps_freq.txt"),
                               sep="\t",header=1,row.names=1)
        } else if (target=="presabs") {
            snps <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_presabs.txt"),
                         sep="\t", row.names=1, header=1)
        } else if (target=="copynum") {
            snps <- read.table(paste0(midas_merge_dir,"/",sp,"/genes_copynum.txt"),
                         sep="\t", row.names=1, header=1)
        } else {
            stop("please specify one of snps, presabs or copynum")
        }

        d <- dist(t(snps), method=distMethod)
        gr <- NULL
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