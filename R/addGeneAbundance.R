#' addGeneAbundance
#' 
#' Add the specified gene copy number to metadata
#' 
#' @param stana stana object
#' @param candSp candidate species
#' @param IDs gene IDs to add
#' @param target KO or genes
#' @param how how to combine multiple IDs
#' @param newCol new column name
#' @param discNumeric convert discrete value to numeric
#' @param disc discretize the abundance by the threshold. function for calculating
#' threshold, like median
#' @param convert conversion such as log10
#' @export
#' @return stana object
addGeneAbundance <- function(stana, candSp, IDs,
    target="KO", how=sum, newCol="gene",
    disc=NULL, discNumeric=TRUE, convert=NULL) {
    if (target=="KO") {
        ints <- intersect(row.names(stana@kos[[candSp]]), IDs)
        subMat <- stana@kos[[candSp]][ints, ]
    } else {
        ints <- intersect(row.names(stana@genes[[candSp]]), IDs)
        subMat <- stana@genes[[candSp]][ints, ]
    }
    if (length(IDs)>1) {
        adda <- apply(subMat, 2, how)    
    } else {
        adda <- subMat
    }
    nm <- names(adda)
    if (!is.null(disc)) {
        thresh <- do.call(disc, list("x"=adda))
        if (discNumeric) {
            adda <- as.numeric(adda > thresh)        
        } else {
            adda <- adda > thresh
        }
        names(adda) <- nm
    }
    meta <- stana@meta
    if (!is.null(convert)) {
        adda <- do.call(convert, list(x=adda))
    }
    meta[[newCol]] <- adda[row.names(stana@meta)]
    stana@meta <- meta
    return(stana)
}

