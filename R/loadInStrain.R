
#' loadInStrain
#' 
#' For InStrain, you should provide `compare` output
#' It can import outputs with --bams options implmented from version 1.6.
#' It computes MAF of var_base per position and return the frequency matrix.
#' Note it considers con_base and var_base in info table for calculation.
#' The function assumes the reads mapped to the default database of inStrain.
#' `skip_pool` option, which loads SNV data, is set to FALSE by default.
#' 
#' @param compare_out_dir output directory of compare.
#' need `output` in the directory
#' @param candidate_species candidate species ID, e.g. `GUT_GENOME000024`
#' @param cl named list of grouping
#' @param just_species if TRUE, returns the vector of species IDs
#' (needs pooled SNV data, otherwise use genome wide compraison table)
#' @param fill_na fill NA in resulting data.frame by -1
#' @param skip_pool skip pooled results loading
#' @return stana object
#' @export
loadInStrain <- function(compare_out_dir,
                         candidate_species,
                         cl=NULL,
                         just_species=FALSE,
                         fill_na=TRUE,
                         skip_pool=TRUE) {

    sample_threshold <- 1
    output_list <- list.files(paste0(compare_out_dir,"/output"))
    stana <- new("stana")
    if (!is.null(cl)) {stana@cl <- cl}
    stana@type <- "InStrain"
    stana@mergeDir <- compare_out_dir
    snps <- list()
    snpsInfoList <- list()

    ctpath <- paste0(compare_out_dir,"/output/",output_list[grepl("comparisonsTable.tsv", output_list)])
    ct <- data.table::fread(ctpath)
    stana@comparisonTable <- ct
    gcpath <- paste0(compare_out_dir,"/output/",output_list[grepl("genomeWide_compare.tsv", output_list)])
    gc <- data.table::fread(gcpath)
    stana@genomeWideCompare <- gc
    scpath <- paste0(compare_out_dir,"/output/",output_list[grepl("strain_clusters.tsv", output_list)])
    sc <- data.table::fread(scpath)
    stana@strainClusters <- sc
    
    if (just_species) {if (skip_pool) {
      return(gc$genome |> strsplit("\\.") |> vapply("[", 1, FUN.VALUE="character") |> unique())
    }}
    stana@ids <- candidate_species

    if (skip_pool) {return(stana)}

    if (sum(grepl("pooled_SNV_data",output_list))!=0) {
        if (!just_species) {
          qqcat("Loading allele count table\n")
          path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data.tsv.gz",output_list)])
          qqcat("Loading the large table...\n")
          tbl <- data.table::fread(path)
          # qqcat("Filtering all-zero positions...\n")
          # tbl <- tbl[apply(tbl, 1, function(x) sum(x[4:7] |> as.numeric())!=0), ]         
        }

        keys_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_data_keys.tsv",output_list)])
        qqcat("Loading key table\n")
        keys <- data.table::fread(keys_path)
        info_path <- paste0(compare_out_dir,"/output/",output_list[grepl("pooled_SNV_info.tsv.gz",output_list)])
        qqcat("Loading info table\n")
        info <- data.table::fread(info_path)

        info$scaffold_position <- paste0(info$scaffold,"|",info$position)
        info$scaffold_position_convar <- paste0(info$scaffold,"|",info$position,"|",info$con_base,">",info$var_base)
        
        varBase <- info$var_base
        names(varBase) <- info$scaffold_position
        conBase <- info$con_base
        names(conBase) <- info$scaffold_position
        
        ## Assuming default database
        sps <- unique(paste0(sapply(strsplit(keys$scaffold, "_"),"[",1),"_",sapply(strsplit(keys$scaffold, "_"),"[",2)))
        if (just_species) {
            return(sps)
        }

        snpsInfoList[[candidate_species]] <- info[grepl(candidate_species, info$scaffold),]
        qqcat("Candidate species: @{candidate_species}\n")
        candKeys <- keys[grepl(candidate_species,keys$scaffold),]$key
        qqcat("  Candidate key numbers: @{length(candKeys)}\n")
        ret <- function(x) {
            smp <- x[1]
            scaff <- x[2]
            scp <- paste0(scaff,"|",as.numeric(x[3]))
            det <- c(smp, scp)
            cb <- x[8]
            vb <- x[9]
            
            read_num <- as.numeric(x[4:7])
            names(read_num) <- names(x[4:7])
            ## Return -1 in zero depth
            if (sum(read_num)==0) {val <- -1} else {
                val <- as.numeric((read_num / sum(read_num))[vb])
            }
            det <- c(det, cb, vb, val)
            det
        }
        mafs <- NULL
        
        sdt <- tbl[tbl$scaffold %in% candKeys,]
        qqcat("  Dimension of pooled SNV table for species: @{dim(sdt)[1]}\n")
        # Only bi-allelic position
        # sdt <- sdt[apply(sdt[,4:7],1,function(x)sum(x==0))==2,]

        nm_sample <- keys$sample
        names(nm_sample) <- keys$key
        nm_sample <- nm_sample[nm_sample!=""]
        nm_scaffold <- keys$scaffold
        names(nm_scaffold) <- keys$key
        nm_scaffold <- nm_scaffold[nm_scaffold!=""]
        
        sdt$sample <- nm_sample[as.character(sdt$sample)]
        sdt$scaffold <- nm_scaffold[as.character(sdt$scaffold)]
        sdt$conBase <- conBase[paste0(sdt$scaffold,"|",sdt$position)]
        sdt$varBase <- varBase[paste0(sdt$scaffold,"|",sdt$position)]
        

        qqcat("Calculating MAF\n")
        maf <- apply(sdt, 1, function(x) ret(x))
        maf <- data.frame(t(maf))
        
        maf <- maf |> `colnames<-`(c("id","scaffold_position","major_allele","minor_allele","maf"))
        maf$maf <- as.numeric(maf$maf)
        maf$scaffold_position <- paste0(maf$scaffold_position,"|",maf$major_allele,">",maf$minor_allele)

        sample_snv <- tidyr::pivot_wider(maf, id_cols=id,
                                 names_from = scaffold_position,
                                 values_from = maf)
        thresh <- apply(sample_snv[2:ncol(sample_snv)], 2, function(x) sum(!is.na(x))) > sample_threshold
        snvDf <- data.frame(sample_snv[,c("id",names(thresh[thresh]))], check.names=FALSE)
        row.names(snvDf) <- snvDf$id
        snvDf$id <- NULL
        snvDf <- data.frame(t(snvDf), check.names=FALSE)
        ## Return the data.frame
        if (fill_na) {
          snvDf[is.na(snvDf)] <- -1
        }
        snps[[candidate_species]] <- snvDf
        stana@snps <- snps
        stana@snpsInfo <- snpsInfoList
        ## Later usage
        stana@snpsInfo[[candidate_species]]$major_allele <- stana@snpsInfo[[candidate_species]]$con_base
        stana@snpsInfo[[candidate_species]]$minor_allele <- stana@snpsInfo[[candidate_species]]$var_base

        stana <- initializeStana(stana,cl)
    } else {
        qqcat("No pooled results present\n")
    }
    return(stana)
}