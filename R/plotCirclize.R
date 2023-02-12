#' plotCirclize
#' 
#' Using `circlize` package, make a circular visualization of 
#' the relationship between genes and corresponding SNVs in genes.
#' @param stana stana object
#' @param genome_id genome ID to plot
#' @param thresh_snp_gene include only genes with at least this number of SNVs
#' @param type one of `point`, `bar`, or `rect` each corresponding to circos functions
#' @param featList features to plot
#' @param cols color of each feature
#' @param featThresh for numeric features, the values above this value will be colored , default to zero.
#' @param cex point size
#' @param textCex text size
#' @param bar_width on circos.barplot
#' @export
#' 
plotCirclize <- function(stana, genome_id, thresh_snp_gene=5, type="point", featList=list(),
                         cols=c("steelblue","tomato","gold","seagreen"),
                         featThresh=0, cex=0.3, textCex=0.5,bar_width=10) {
    qqcat("Type is @{stana@type}\n")
    info <- stana@snpsInfo[[candSp]]

    if (length(featList)==0) {
        qqcat("Features not provided, default to sample_counts\n")
        if (stana@type=="MIDAS2") {
          sample_counts <- info$sample_counts
        } else if (stana@type=="InStrain") {
          sample_counts <- info$sample_detections # Probably 5x_detections
        }
        names(sample_counts) <- row.names(info)
        featList[["sample_counts"]] <- sample_counts
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
      qqcat("  @{i}: @{min(tmp_info$position)} - @{max(tmp_info$position)}, number of position: @{nrow(tmp_info)}\n")
    }
    
    qqcat("Genome ID: @{genome_id}\n")

    ## Deleting "None" and gene ids with only one position
    circ_plot <- subset(info, info$genome_id==genome_id)
    for (nm in names(featList)) {
        circ_plot[[nm]] <- featList[[nm]][row.names(circ_plot)]
        circ_plot <- circ_plot[!is.na(circ_plot[[nm]]),]
    }

    ## Sector will be gene_id
    circ_plot <- subset(circ_plot, gene_id!="None")
    gene_num <- table(circ_plot$gene_id)
    circ_plot <- subset(circ_plot, !gene_id %in% names(gene_num[gene_num<thresh_snp_gene]))
    
    qqcat("Included position: @{dim(circ_plot)[1]}\n")


    circos.clear() ## If have one
    circos.par(cell.padding=c(0.02,0,0.02,0), gap.degree=0.5)
    circos.initialize(factors = circ_plot$gene_id,
                      x=as.numeric(circ_plot$position))
    for (i in seq_along(names(featList))) {
      
        nm <- names(featList)[i]
        colfunc<-colorRampPalette(c("blue","red"))
        circ_plot$interval <- findInterval(circ_plot[[nm]], sort(circ_plot[[nm]]))
        # circ_plot$sq_position <- 1:nrow(circ_plot)
        circ_plot$color <- colfunc(nrow(circ_plot))[circ_plot$interval]
        
        if (i==1) {
          ## Plot text bending inside
          if (type=="point") {
            circos.track(
              factors=circ_plot$gene_id,
              track.margin = c(0, mm_h(2)),
              x=as.numeric(circ_plot$position),
              y=as.numeric(circ_plot[[nm]]),
              ylim=c(-1.1,1.1),
              panel.fun = function(x, y) {
                circos.text(CELL_META$xcenter,
                            CELL_META$cell.ylim[2] + mm_y(2),
                            CELL_META$sector.index,
                            facing="bending.inside",
                            cex=textCex)
                circos.points(x, y, pch=19, cex=cex,
                              # col = col)
                              col = ifelse(y > featThresh, cols[i], "black"))
                
                
              })
          } else if (type=="rect") {
            
            circos.track(ylim = c(0, 1),
                         factors=circ_plot$gene_id,
                         x=circ_plot$position,
                         y=circ_plot$color,
                         panel.fun = function(x, y) {
                           
                           circos.text(CELL_META$xcenter,
                                       CELL_META$cell.ylim[2] + mm_y(2),
                                       CELL_META$sector.index,
                                       facing="bending.inside",
                                       cex=textCex)
                           
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
                
                circos.text(CELL_META$xcenter,
                            CELL_META$cell.ylim[2] + mm_y(2),
                            CELL_META$sector.index,
                            facing="bending.inside",
                            cex=textCex)
                
                circos.barplot(value=y, pos=x, bar_width = bar_width,
                               border=ifelse(y > featThresh, cols[i], "black"),
                               col = ifelse(y > featThresh, cols[i], "black"))
              })            
          }
        } else {
          
          if (type=="point") {
              circos.track(
                       factors=circ_plot$gene_id,
                       track.margin = c(0, mm_h(2)),
                       x=as.numeric(circ_plot$position),
                       y=as.numeric(circ_plot[[nm]]),
                       ylim=c(-1.1,1.1),
                       panel.fun = function(x, y) {
                         circos.points(x, y, pch=19, cex=cex,
                                         # col = col)
                                         col = ifelse(y > featThresh, cols[i], "black"))
                           
                           
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
                              print(x[-n])
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
                                     border=ifelse(y > featThresh, cols[i], "black"),
                                    col = ifelse(y > featThresh, cols[i], "black"))
                  })            
          }
        }
    }
    circos.clear()
}
