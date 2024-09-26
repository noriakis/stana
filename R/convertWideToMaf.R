#' convertWideToMaf
#' 
#' should be six column layout with the column corresponding to
#' "species_id","snv_id","ref_allele","alt_allele","ref_count","alt_count"
#' @param df data frame of six columns
#' @return species-named list of MAF data.frame
#' @export
convertWideToMaf <- function(df) {
    colnames(df) <- c("species_id","snv_id","ref_allele","alt_allele","ref_count","alt_count")
    df$ref_count <- as.numeric(df$ref_count)
    df$alt_count <- as.numeric(df$alt_count)

    df$count <- df$ref_count + df$alt_count
    df$maf <- apply(df, 1, function(x) {
        if (as.numeric(x[5]) > as.numeric(x[6])) {
            as.numeric(x[6]) / as.numeric(x[7])
        } else {
            as.numeric(x[5]) / as.numeric(x[7])
        }
    })
    df$major <- apply(df, 1, function(x) {
        if (as.numeric(x[5]) > as.numeric(x[6])) {
            x[3]
        } else {
            x[4]
        }
    })
    df$minor <- apply(df, 1, function(x) {
        if (as.numeric(x[5]) > as.numeric(x[6])) {
            x[4]
        } else {
            x[3]
        }
    })
    sp <- unique(df$species_id)
    ret <- lapply(sp, function(s) {
        tmp <- df %>% dplyr::filter(.data$species_id == s)
        tmp <- tmp[,c("snv_id", "maf")]
        row.names(tmp) <- tmp$snv_id
        tmp$snv_id <- NULL
        tmp
    })
    names(ret) <- sp

    summary <- lapply(sp, function(s) {
        tmp <- df %>% dplyr::filter(.data$species_id == s)
        tmp <- tmp[,c("snv_id", "major", "minor")]
        row.names(tmp) <- tmp$snv_id
        tmp$snv_id <- NULL
        colnames(tmp) <- c("major_allele","minor_allele")
        tmp$species_id <- s
        tmp
    })
    names(ret) <- sp
    names(summary) <- sp

    return(list(ret, summary))

}

#' combineMaf
#' combine the maf table produced by convertWideToMaf()
#' @param mafs named list of the results of convertWideToMaf()
#' @export
#' @return stana object
combineMaf <- function(mafs) {
  if (is.null(names(mafs))) {stop("The list must be named")}
  sps <- do.call(union, lapply(names(mafs), function(s) {
    ## Species
    names(mafs[[s]][[1]])
  }))
  snvls <- lapply(sps, function(s) {
    tmpdf <- do.call(rbind, lapply(names(mafs), function(sample) {
      tmp <- mafs[[sample]][[1]][[s]]
      if (is.null(tmp)) {return(tmp)}
      tmp$sample <- sample
      tmp$snv_id <- row.names(tmp)
      return(tmp)
    }))
    tbl <- tidyr::pivot_wider(tmpdf,
                       id_cols="snv_id",
                       names_from="sample",
                       values_from = "maf")
    tbl <- data.frame(tbl)
    row.names(tbl) <- tbl$snv_id
    tbl$snv_id <- NULL
    tbl
  })
  
  snvsu <- lapply(sps, function(s) {
    tmpdf <- do.call(rbind, lapply(names(mafs), function(sample) {
      tmp <- mafs[[sample]][[2]][[s]]
      if (is.null(tmp)) {return(tmp)}
      tmp$sample_name <- sample
      return(tmp)
    }))
    tbl <- data.frame(tmpdf)
    tbl
  })
  
  names(snvls) <- sps
  stana <- new("stana")
  stana@type <- "manual"
  stana@snps <- snvls
  stana@snpsInfo <- snvsu
  stana
}