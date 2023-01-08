#'
#' doAdonis
#' 
#' Perform PERMANOVA on distance matrix based 
#' on SNV frequency or gene matrix.
#' Named list of samples are to be provided.
#' Note that this function currently performs 
#' comparison using one category.
#' for complex design, one should do it manually.
#' 
#' @param stana stana object
#' @param specs species to be examined
#' @param cl named list of samples
#' @param target snps, presabs, or copynum
#' @param formula pass to adonis2, specify distance matrix as d.
#' @param distMethod distance method passed to dist() (default, manhattan)
#' @param maj major allele distance
#' @param ... parameters passed to adonis2
#' @importFrom vegan adonis2
#' @importFrom stats as.formula dist
#' @importFrom utils read.table
#' @export
doAdonis <- function(stana, specs, cl,
    target="snps", formula=NULL,
    distMethod="manhattan",
    maj=FALSE,
    ...) {
      for (sp in specs){
        qqcat("Performing adonis in @{sp}\n")
        if (target=="snps") {
            snps <- stana@snps[[sp]]
        } else {
            snps <- stana@genes[[sp]]
        }
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
        stana@adonisList[[sp]] <- adores
      }
      stana
}