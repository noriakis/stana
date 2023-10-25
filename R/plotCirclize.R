#' plotCirclize
#' 
#' Using `circlize` package, make a circular visualization of 
#' the relationship between genes and corresponding SNVs in genes.
#' @param stana stana object
#' @param candSp candidate species ID
#' @param genomeId genome ID to plot
#' @param thresh_snp_gene include only genes with at least this number of SNVs
#' @param featList features to plot
#' @param cols color of each feature
#' @param controlColor control (not reaching `featThresh`) color in points and barplot
#' @param featThresh for numeric features, the values above this value will be colored , default to zero.
#' @param cex point size
#' @param MAF plot MAF for each group
#' @param cl cluster information (named list)
#' @param textCex text size
#' @param showGeneName show gene ID on the plot
#' @param bar_width on circos.barplot
#' @param contPalette palette for continuous scale
#' @param discPalette palette for discrete scale
#' @param returnRawDf return the raw data frame used in the plot
#' @param featCircos which type of circos to use in the corresponding feature
#' @param include_gene include these gene IDs, ignoring the thesh_snp_gene argument
#' @param only_genome_id show only the genome IDs
#' @importFrom circlize CELL_META circos.clear circos.par circos.initialize circos.rect circos.barplot circos.track circos.text circos.points
#' @return draw circlize plot
#' @export
#' 
plotCirclize <- function(stana, candSp, genomeId, include_gene=NULL,
                         thresh_snp_gene=5, featList=list(),
                         cols=c("tomato","steelblue","gold","seagreen"),
                         controlColor="steelblue", MAF=FALSE, cl=NULL,
                         featCircos=list(), showGeneName=TRUE,
                         returnRawDf=FALSE, only_genome_id=FALSE,
                         featThresh=0, cex=0.3, textCex=0.5,bar_width=10,
                         contPalette=c("steelblue",
                          "tomato"), discPalette="Dark2") {
    qqcat("Type is @{stana@type}\n")
    if (is.null(stana@snpsInfo[[candSp]])) {stop("No SNV info available")}
    info <- stana@snpsInfo[[candSp]]
    if (MAF) {
      maf <- stana@snps[[candSp]]
      if (is.null(cl)) {cl <- stana@cl}
      ## Subset to MAF row.names
      info <- info[row.names(maf),]
      for (clnm in names(cl)) {
        inc <- intersect(colnames(maf),cl[[clnm]])
        tmpMaf <- maf[, inc]
        ## For considering MIDAS2, -1 will be deleted
        featList[[clnm]] <- apply(tmpMaf, 1, function(x) mean(x[x != -1]))
        featCircos[[clnm]] <- "rect"
      }
    }

    if (length(featList)==0) {
        qqcat("Features not provided, default to sample_counts\n")
        ## Choosing default type for each pipeline
        if (stana@type=="MIDAS2") {
          sample_counts <- info$sample_counts
        } else if (stana@type=="InStrain") {
          sample_counts <- info$sample_detections # Probably 5x_detections
        }
        names(sample_counts) <- row.names(info)
        featList[["sample_counts"]] <- sample_counts
        featCircos[["sample_counts"]] <- "rect"
    } else {
      qqcat("Number of features: @{length(featList)}\n")
    }

    if (stana@type=="MIDAS2") {
      info$position <- as.numeric(sapply(strsplit(row.names(info),"\\|"),"[",4))
      info$genome_id <- sapply(strsplit(row.names(info),"\\|"),"[",3)
    } else if (stana@type=="InStrain") {
      info$genome_id <- info$scaffold
    }

    qqcat("Genome ID in SNV information:\n")
    for (i in unique(info$genome_id)) {
      tmp_info <- subset(info, genome_id==i)
      qqcat("  @{i}: @{min(tmp_info$position)} - @{max(tmp_info$position)}",
      	", number of position: @{nrow(tmp_info)}\n")
    }
    
    if (only_genome_id) {return(1)}
    
    qqcat("Genome ID: @{genomeId}\n")

    ## Deleting "None" and gene ids with only one position
    info$genome_id <- info$genome_id |> strsplit("_") |> 
                                        vapply("[", 1, FUN.VALUE="character")
    circ_plot <- subset(info, info$genome_id==genomeId)
    if (length(featList)>0) {
      for (nm in names(featList)) {
          circ_plot[[nm]] <- featList[[nm]][row.names(circ_plot)]
          circ_plot <- circ_plot[!is.na(circ_plot[[nm]]),]
      }      
    }
    ## Sector will be gene_id
    circ_plot <- subset(circ_plot, circ_plot$gene_id!="None")
    gene_num <- table(circ_plot$gene_id)
    if (is.null(include_gene)) {
      circ_plot <- subset(circ_plot, !circ_plot$gene_id %in% names(gene_num[gene_num<thresh_snp_gene]))
    } else {
      circ_plot <- subset(circ_plot, circ_plot$gene_id %in% include_gene)
    }
    qqcat("Included position: @{dim(circ_plot)[1]}\n")
    if (dim(circ_plot)[1]==0) {stop("No position available, finishing ...")}
    if (returnRawDf) {return(circ_plot)}
    circos.clear() ## If have one
    circos.par(cell.padding=c(0.02,0,0.02,0),
      gap.degree=0.5)
    circos.initialize(factors = circ_plot$gene_id,
                      x=as.numeric(circ_plot$position))
    for (i in seq_along(names(featList))) {
        nm <- names(featList)[i]
        type <- featCircos[[nm]]
        if (is.numeric(featList[[nm]])) {
          colfunc<-colorRampPalette(c(contPalette[1],contPalette[2]))
          circ_plot$interval <- findInterval(circ_plot[[nm]], sort(circ_plot[[nm]]))
          # circ_plot$sq_position <- 1:nrow(circ_plot)
          circ_plot$color <- colfunc(nrow(circ_plot))[circ_plot$interval]
        } else {
          colfunc<-colorRampPalette(brewer.pal(length(unique(circ_plot[[nm]])), discPalette))
          discColors <- colfunc(length(unique(circ_plot[[nm]])))
          names(discColors) <- unique(circ_plot[[nm]])
          circ_plot$color <- sapply(circ_plot[[nm]], function(x) discColors[x])
        }


        if (i==1) {
          ## Plot text bending inside if i==1
          if (type=="point") {
            circos.track(
              factors=circ_plot$gene_id,
              track.margin = c(0, mm_h(2)),
              x=as.numeric(circ_plot$position),
              y=as.numeric(circ_plot[[nm]]),
              panel.fun = function(x, y) {
                if (showGeneName){
                  circos.text(CELL_META$xcenter,
                              CELL_META$cell.ylim[2] + mm_y(2),
                              CELL_META$sector.index,
                              facing="bending.inside",
                              cex=textCex)                  
                }
                circos.points(x, y, pch=19, cex=cex,
                              # col = col)
                              col = ifelse(y > featThresh, cols[i], controlColor))
                
                
              })
          } else if (type=="rect") {
            
            circos.track(ylim = c(0, 1),
                         factors=circ_plot$gene_id,
                         x=circ_plot$position,
                         y=circ_plot$color,
                         panel.fun = function(x, y) {
                           if (showGeneName) {
                             circos.text(CELL_META$xcenter,
                                         CELL_META$cell.ylim[2] + mm_y(2),
                                         CELL_META$sector.index,
                                         facing="bending.inside",
                                         cex=textCex)                            
                           }                           
                           xlim = CELL_META$xlim
                           ylim = CELL_META$ylim
                           n = length(x)
                           ## Position are squeezed, unlike point
                           circos.rect(x[-n], rep(ylim[1], n-1),
                                       x[-1], rep(ylim[2], n-1),
                                       col = y,
                                       border = NA)
                         })
          } else {
            circos.track(
              factors=circ_plot$gene_id,
              track.margin = c(0, mm_h(2)),
              x=as.numeric(circ_plot$position),
              y=as.numeric(circ_plot[[nm]]),
              panel.fun = function(x, y) {
                if (showGeneName) {
                  circos.text(CELL_META$xcenter,
                              CELL_META$cell.ylim[2] + mm_y(2),
                              CELL_META$sector.index,
                              facing="bending.inside",
                              cex=textCex)                  
                }                
                circos.barplot(value=y, pos=x, bar_width = bar_width,
                               border=ifelse(y > featThresh, cols[i], controlColor),
                               col = ifelse(y > featThresh, cols[i], controlColor))
              })            
          }
        } else {
          
          if (type=="point") {
              circos.track(
                       factors=circ_plot$gene_id,
                       track.margin = c(0, mm_h(2)),
                       x=as.numeric(circ_plot$position),
                       y=as.numeric(circ_plot[[nm]]),
                       # ylim=c(-1.1,1.1),
                       panel.fun = function(x, y) {
                         circos.points(x, y, pch=19, cex=cex,
                                         # col = col)
                                         col = ifelse(y > featThresh, cols[i], controlColor))
                           
                           
                       })
          } else if (type=="rect") {
  
              circos.track(ylim = c(0, 1),
                           factors=circ_plot$gene_id,
                           x=circ_plot$position,
                           y=circ_plot$color,
                           panel.fun = function(x, y) {
                              xlim = CELL_META$xlim
                              ylim = CELL_META$ylim
                              n = length(x)
                              ## Position are squeezed, unlike point
                              circos.rect(x[-n], rep(ylim[1], n-1),
                                          x[-1], rep(ylim[2], n-1),
                                          col = y,
                                          border = NA)
              })
          } else {
              circos.track(
                  factors=circ_plot$gene_id,
                  track.margin = c(0, mm_h(2)),
                  x=as.numeric(circ_plot$position),
                  y=as.numeric(circ_plot[[nm]]),
                  panel.fun = function(x, y) {
                      circos.barplot(value=y, pos=x, bar_width = bar_width,
                                     border=ifelse(y > featThresh, cols[i], controlColor),
                                    col = ifelse(y > featThresh, cols[i], controlColor))
                  })            
          }
        }
    }
    circos.clear()
}
