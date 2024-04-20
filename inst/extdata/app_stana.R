## Interactive application version 0.99.0
## Dependencies
library(stana)
## Shiny-related
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(waiter)
library(dplyr)
library(ggtree)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggraph)
library(ggkegg)
library(ggnewscale)
library(ggtreeExtra)
library(ggstar)
library(MKmisc)
library(purrr)
library(castor)
library(pillar)
library(ComplexHeatmap)

## Random palette
pal <- sample(RColorBrewer::brewer.pal.info %>% row.names(), 1)

set.seed(1)

cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))

## NNLM preferred (fast)
if (requireNamespace("NNLM")) {
	library(NNLM)
	nnlm <- TRUE
} else {
	if (requireNamespace("NMF")) {
		library(NMF)
		nnlm <- FALSE
	} else {
		stop("NNLM or NMF package should be installed.")
	}
}


## Load the necessary data
if (!file.exists("ko_pathway.rda")) {
    kop <- data.table::fread("https://rest.kegg.jp/link/ko/pathway", header=FALSE)
    save(file="ko_pathway.rda", kop, compress="xz")
} else {
    load(file="ko_pathway.rda")
}
if (!file.exists("pathway_name.rda")) {
    nmls <- data.table::fread("https://rest.kegg.jp/list/pathway/ko", header=FALSE)
    save(file="pathway_name.rda", nmls, compress="xz")
} else {
    load(file="pathway_name.rda")
}

list_of_species <- NULL

## Function to reset the results
reset <- function(output) {
    output$showKEGGPlot <- NULL
    output$path_selector <- NULL
    output$species_name <- NULL
    output$nSample <- NULL
    output$adonis <- NULL
    output$showPlot <- NULL
    output$showSigTree <- NULL
}
reset_enr <- function(output) {
    output$showEnr <- NULL
    output$showGenePlot <- NULL
}

## Moderated T statistics calculation based on MKmisc implementation
mod.T <- function(mat, l1, l2) {
    if (dim(mat)[1]==0) {return(NULL)}
    ordered.mat <- mat[, c( intersect(colnames(mat), l1), intersect(colnames(mat), l2) )]
    gr <- c( rep("l1", length(intersect(colnames(mat), l1))), rep("l2", length(intersect(colnames(mat), l2))) )
    res <- MKmisc::mod.t.test(as.matrix(ordered.mat), gr)
    modt <- res$t
    names(modt) <- row.names(res)
    return(modt)
}

## Return Function x Species heatmap
return_all_gsea_heatmap <- function(list_of_gsea, sel="enrichmentScore",
                                    retTable=FALSE, title="GSEAMat") {
    matl <- lapply(names(list_of_gsea), function(x) {
        list_of_gsea[[x]] |> dplyr::select(pathway, sel)
    })
    hm <- matl |>
        purrr::reduce(full_join, by = "pathway") |>
        `colnames<-`(c("pathway",names(list_of_gsea))) |>
        data.frame(check.names=FALSE)

    row.names(hm) <- hm$pathway
    hm$pathway <- NULL
    hm2 <- hm |> mutate_all(~as.numeric(.))
    changer <- nmls$V2 |> setNames(nmls$V1)
    row.names(hm2) <- changer[row.names(hm2)]
    if (retTable) {return(hm2)}
    Heatmap(hm2 |> as.matrix(), cluster_rows = FALSE,
            column_title=title,
            column_names_max_height = unit(10, "cm"))        
}

## Dataset definition
dats <- list.files("data")

## Stana object is stored in RDA object
datasets <- dats[grepl(pattern = ".rda", dats)] %>% strsplit(".rda") %>% vapply("[", 1, FUN.VALUE="a")
names(datasets) <- datasets

first_choice <- NA

ui <-  fluidPage(
    useWaiter(),
    id="main",
    useShinyjs(),
    title = "stana shiny interactive application",
    tabsetPanel(id = "main_tab",
        tabPanel("Tree",
            fluidRow(
                column(12, align="center",
                div(plotOutput("showPlot", height="100%")),
                 ## Adonis output
                div(DT::dataTableOutput("adonis"),style= "font-size: 75%"),
                ),
            hr(),
            fluidRow(
                ## Species selection column
                column(3, align="center",
                    selectInput("site", label = "Dataset", choices = datasets),
                    actionButton("inspect", "Order by the number of samples", class = "btn-success"),
                    p("Output number of samples in the trees in Inspect tab"),
                    shiny::hr(),
                    selectInput("species","Select species", choices = list_of_species),
                    numericInput("all","Multiple",1),
                    numericInput("colnum", "Number of column",3)
                ),

                column(3, align="center",
                    selectInput('cv', "Clinical variables",
                        choices = NULL, multiple=TRUE),
                    uiOutput("col")
                ),

                column(3, align="center",
                    h4("Visualization options"),
                    numericInput("width","width",value=1000),
                    numericInput("height","height",value=1000),
                    numericInput("size","Size",value=3),
                    selectInput("layout","Layouts",
                    choices = c("rectangular","fan","slanted", "circular")),
                    switchInput("cladogram", label="Cladogram"),
                    switchInput("use_geom_fruit", label="Use geom_fruit"),
                    switchInput("use_ko_dist", label="Use KO"),
                    numericInput("rank", "Rank", value=2),
                    switchInput("use_nmf", label="Use NMF"),
                    switchInput("use_pairwise", label="Pairwise PERMANOVA"),
                    actionButton("show", "Show", class = "btn-success"),
                    hr(),
                    br(),
                    numericInput("thresh", "Sample filtering threshold (fraction of whole metadata,
                        if > 1, use the thresholding number)", value=0.2),
                    br(),
                    selectInput("hmcell", "Heatmap cell", choices=c("pseudoF","P","R")),
                    br(),
                    actionButton("alladonis", "Perform PERMANOVA for all", class = "btn-success"),
                    hr(),
                    numericInput("save_width","Save width",value=10),
                    numericInput("save_height","Save height",value=10),
                    downloadButton("downloadPlot", "Download plot"),
                    downloadButton("downloadHeatmap", "Download heatmap"),
                    downloadButton("downloadMat", "Download matrix")
                ),

                column(3, align="center",
                    h4("Gene options"),
                    numericInput("kmeans", "K-means", value=10),
                    numericInput("reduce", "Reduce", value=100),
                    switchInput("lfc", label="T-statistics (must be two levels)"),
                    p("If not T-stat, sum values of KO abundance across sample",
                      "is passed to EA and KEGG"),
                    switchInput("gsea", label="Perform GSEA instead of ORA"),
                    p("ORA options"),
                    switchInput("pvalAdjustMethod", label="Adjust method", onLabel="BH", offLabel="none"),
                    switchInput("updown", label="Up or Down when ORA", onLabel="Up", offLabel="Down"),
                    actionButton("genes", "Check CNV",
                        class = "btn-success"),
                    hr(),
                    p("GSEA / Species"),
                    actionButton("allgsea", "Produce GSEA matrix (experimental)", class="btn-success"),
                    downloadButton("downloadGSEAMat", "Download GSEA matrix")
                )
            )
        ),
        hr(),
        fluidRow(
            column(12, align="center",
                verbatimTextOutput("log", placeholder = TRUE)
            )
        )
    ), ## End of Main Panel
    tabPanel("CNVs",
        div(plotOutput("showGenePlot", height = "100%"),
            style="overflow-y: scroll;overflow-x: scroll;"
        ),
        conditionalPanel(
            condition="input.genes>=1",
            h4("Enrichment analysis results"),
            div(DT::dataTableOutput("showEnr"),
            style="height:500px; overflow-y: scroll;overflow-x: scroll;")
        ),
        conditionalPanel(
            condition="input.genes>=1",
            selectInput("path_selector", "Pathway", choices="")
        ),
        conditionalPanel(
            condition="input.genes>=1",
            switchInput("ggraph", label="Use ggraph (only on global map)"),
            actionButton("kegg", "Check KEGG", class = "btn-success")
        )
    ), ## End of CNV panel
    tabPanel("KEGG PATHWAY",
        div(plotOutput("showKEGGPlot", height = "100%")),
        div(DT::dataTableOutput("KOSTAT"),style= "font-size: 75%"),
        downloadButton("downloadKEGGPathwayPlot", "Download pathway plot")
    ),
    tabPanel("Inspect",
        column(12, align="center",
            div(plotOutput("showSigTree", height = "100%")),
            div(DT::dataTableOutput("stats"),style="overflow-y: scroll;overflow-x: scroll;"),
            downloadButton("downloadStatTSV", "Download TSV")
        )
    ),
    tabPanel("Search",
        h1("Search species"),
        textInput("searchsp",label="Species name (e.g., Veillonella dispar)"),
        actionButton("search", "Search", class = "btn-success"),
        DT::dataTableOutput("search"),
        hr(),
        selectInput("search_input", "List", ""),
        actionButton("search_each", "Show species tree",
        class = "btn-success")
    )
))

server <- function(input, output, session) {
    values <- reactiveValues()
    ## Add metadata
    col_names <- reactive(input$cv)
    output$col <- renderUI({
        map(col_names(), ~ switchInput(.x,
            label=.x, onLabel = "numeric",
            offLabel="factor",
            value = isolate(input[[.x]])))
    })

    ## Logging functions
    a <- reactiveValues(logOutput="[log]")
    output$log <- renderText({ a$logOutput })
    updateLog <- function(text){
        a$logOutput <- paste(a$logOutput, paste0(Sys.time(), ": ", text), sep = "\n")
    }

    ## When changing the dataset, load metadata and relevant species
    observeEvent(input$site, {
        reset(output)
        load(paste0("data/",input$site,".rda"))
        values$stana <- stana
        all_id <- stana@ids
        if (length(all_id)==0) {
        	showNotification("No species available in the dataset")
        	return(1)
        }
        values$all_id <- all_id
        values$all_names <- stana@names
        values$rev_names <- names(stana@names) %>% setNames(stana@names)

        ## This will be shown to users
        list_of_species <- as.character(stana@names)

        updateLog(paste0("Changing site to ", input$site))
        cat_subtle("# Changing site to ", input$site, "\n", sep="")

        updateSelectInput(session, "species", choices=list_of_species)
        values$current_list_of_species <- list_of_species

        if (input$search_input!="") {
            search_for_sp <- values$current_list_of_species
            update_for_sp <- search_for_sp[grepl(values$search_spid, search_for_sp)]
            updateSelectInput(session, "species", selected=update_for_sp)
            updateSelectInput(session, "search_input", selected="")
        }
        ## Metadata loading
        if (dim(stana@meta)[1]==0) {
            shiny::showNotification("No metadata found")
            return(1)
        }
        ## Here, we just update column name to select input,
        ## so no mod on metadata
        meta <- stana@meta
        ## Load metadata onto reactive
        values$meta <- meta

        updateSelectInput(session, 'cv',
            choices = colnames(meta),
            # server = TRUE,
            selected = first_choice[input$site])
    })


    ## Output statistics of the tree and ...
    observeEvent(input$inspect,{
        waiter_show(html = spin_3(), color = "white")
        updateLog("Counting samples ...")
        statsdf <- NULL
        withProgress(
            for (i in values$all_id) {
                tre <- values$stana@treeList[[i]]
                statsdf <- rbind(statsdf,
                                 c(i, values$all_names[i], length(tre$tip.label)))
            },
            message="Reading trees ..."
        )
        statsdf <- statsdf |> data.frame() |> `colnames<-`(c("ID","Species","Number"))
        statsdf$Number <- as.numeric(statsdf$Number)
        updateTabsetPanel(session, "main_tab", selected = "Inspect")
        ordered <- statsdf[order(statsdf$Number, decreasing=TRUE),]
        output$stats <- DT::renderDataTable(ordered)
        values$stats <- ordered
        
        ## Update based on sampled number
        list_of_species_ordered <- ordered$Species
        values$ordered <- ordered$ID
        values$orderedTax <- list_of_species_ordered
        updateSelectInput(session, "species", choices=list_of_species_ordered)
        waiter_hide()}
    )



    ## Search the species across dataset
    observeEvent(input$search, {
        searchlet <- list(); k <- 1;
        if (input$searchsp=="") {return(1);}
        updateLog(paste0("Searching for ", input$searchsp))
        for (i in names(datasets)) {
            searchsp <- input$searchsp
            load(paste0("data/",i,".rda"))
            all_id <- stana@ids
            all_names <- stana@names
            list_of_species <- all_names

            ## Actual search here, grepl
            candidate <- list_of_species[ grepl(searchsp, list_of_species) ]
            if (length(candidate)==0) {} else {
                for (j in candidate) {
                    tmptax <- values$rev_names[j]
                    tre <- stana@treeList[[tmptax]]
                    searchlet[[k]] <- c(i, j, length(tre$tip.label))
                    k <- k + 1
                }
            }
        }
        if (length(searchlet)==0) {return(1);}
        cat_subtle("# Candidate: ", length(candidate), "\n", sep="");
        searchlet <- do.call(rbind, searchlet) |>
            data.frame() |>
            `colnames<-`(c("dataset","taxonomy","number"))
        searchlet$number <- as.numeric(searchlet$number)
        listed <- paste0(searchlet$dataset, ":", searchlet$taxonomy)
        output$search <- DT::renderDataTable(searchlet)
        updateSelectInput(session, "search_input", choices = listed)
    }) ## End of search across dataset

    ## Update the panels based on searched species
    observeEvent(input$search_each, {
        candsp <- input$search_input
        ds <- strsplit(candsp, ":") |> sapply("[",1)
        spid <- strsplit(candsp, ":") |> sapply("[",2)
        datadir <- datasets[[ds]]
        values$search_spid <- spid
        updateTabsetPanel(session, "main_tab", selected = "Tree")
        updateSelectInput(session, "site", selected=datasets[[ds]])

        search_for_sp <- values$current_list_of_species
        update_for_sp <- search_for_sp[grepl(spid, search_for_sp)]
        updateSelectInput(session, "species", selected=update_for_sp)
        updateSelectInput(session, "search_input", selected="")
    })


    ## All adonis: Perform PERMNOVA for all species 
    ## [CAUTION] even when pairwise option is enabled
    ## pairwise calculation will not be performed in this part

    observeEvent(input$alladonis, {
        ## Delete the significant trees and stats tables 
        output$showSigTree <- NULL
        output$stats <- NULL
        waiter_show(html = spin_3(), color = "white")
        if (length(input$cv)>1) {
            updateLog(paste0("Multiple items selected"))
        }
        trels <- list(); metals <- list();
        updateLog(paste0("Performing PERMANOVA multiple times"))
        ps <- list()
        for (sp in values$current_list_of_species) {
            sp_id <- values$rev_names[sp]

            if (input$use_ko_dist) {
                cat_subtle("# Using KO table to calculate distance between samples\n")
                ## Load the precalculated KO (or gene) table
                ko_df_filt <- values$stana@kos[[sp_id]]
                if (length(colnames(ko_df_filt))<3) {
                    cat_subtle("# Too few samples ", sp_id, "\n")
                    next
                }
                dist_mat <- dist(t(ko_df_filt))
                if ("TRUE" %in% as.character(names(table(is.na(dist_mat))))) {
                    next
                }
                tre <- phangorn::upgma(dist_mat)
            } else {
                tre <- values$stana@treeList[[sp_id]]
            }
            
            trels[[sp_id]] <- tre
            labnum <- length(tre$tip.label)

            ## Processing if tree has 0 or 1 sample
            if (is.null(tre$tip.label) |
                length(tre$tip.label)==1) {
                cat_subtle("# ", paste0(sp_id,": Not enough samples are profiled"), "\n", sep="")
            } else {
                meta <- values$stana@meta
                if (input$thresh>=1) {
                    curthre <- input$thresh
                } else {
                    curthre <- dim(meta)[1]*input$thresh
                }
                if ( curthre > labnum) {
                    cat_subtle("# ", paste0(sp_id, ": Filtered by the threshold < ",
                                   dim(meta)[1]*input$thresh), "\n", sep="");
                    next
                }
                ## If tree labels are not in metadata
                ## (which should not happen if exported from stana)
                for (i in tre$tip.label) {
                    if (!i %in% row.names(meta)) {
                        cur <- row.names(meta)
                        meta <- rbind(meta, rep(NA, ncol(meta)))
                        row.names(meta) <- c(cur, i)
                    }
                }

                ## Preprocess
                for (tmp_cv in input$cv) {
                    change_cv <- meta[[tmp_cv]]
                    change_cv[change_cv==""] <- NA
                    change_cv[change_cv=="-"] <- NA
                    change_cv[is.infinite(change_cv)] <- NA

                    if (input[[tmp_cv]]) {
                        meta[[tmp_cv]] <- as.numeric(change_cv)
                    } else  {
                        meta[[tmp_cv]] <- as.factor(change_cv)
                    }
                }
                meta <- meta[tre$tip.label, ]
                if (is.na(meta[tre$tip.label, tmp_cv]) |> sum() == dim(meta)[1]) {
                    cat_subtle("# One of the covariates are all NA\n")
                    updateLog(paste0("One of the covariates are all NA in ", sp," exiting"));
                    next
                }


                if (sum(input$cv %in% colnames(meta))==0) {
                    shiny::showNotification("Selected column is not",
                                            " present in metadata")
                    # waiter_hide();
                    next
                }

                ## Test PERMANOVA
                if (!input$use_ko_dist) {
                    dist_mat <- as.dist(ape::cophenetic.phylo(tre))
                }
                subset_meta <- meta[tre$tip.label, ]

                ## [CAUTION] If NA is present in metadata, they are omitted!
                for (tmp_cv in input$cv) {
                    if ((subset_meta[[tmp_cv]] |> na.omit() |> unique() |> length()) == 1) {
                        updateLog("One of the covariates have one level only")
                        next
                    }
                }

                if (length(input$cv)==1) {
                    adonis_meta <- subset_meta[input$cv] |> as.data.frame()
                } else {
                    adonis_meta <- subset_meta[input$cv]
                }
                
                metals[[sp_id]] <- adonis_meta

                tryCatch(
                    {
                        cat_subtle("Doing ADONIS for ", sp, "\n", sep="")
                        withProgress(
                            adonis_res <- vegan::adonis2(as.formula(paste0("dist_mat ~  .")),
                                                         data=adonis_meta, na.action = na.omit),
                            message=paste0("Doing PERMANOVA on ", sp)
                        )
                    },
                    error=function(e) {print(e); adonis_res <- NULL}
                )
                if (!is.null(adonis_res)) {
                    if (length(input$cv)==1) {
                        ps[[sp_id]] <- c(sp,
                            labnum,
                            adonis_res$R2[1],
                            adonis_res$F[1],
                            adonis_res$`Pr(>F)`[1])
                    } else {
                        ## Save all to list
                        ps[[sp_id]] <- list(sp, labnum,
                                            adonis_res$Df %>% tail(1),
                                            adonis_res$R2[seq_len(length(input$cv))],
                                            adonis_res$F[seq_len(length(input$cv))],
                                            adonis_res$`Pr(>F)`[seq_len(length(input$cv))]) |>
                            setNames(c("sp","numsample","totdf","r2","F","p"))
                    }
                }
            }

        }# SPECIES LOOP

        if (length(ps)!=0) {
            if (length(input$cv)==1) {
                dfs <- do.call(rbind, ps) |>
                    data.frame() |>
                    `colnames<-`(c("species","sample","R2","F","p"))
                dfs <- dfs[!is.na(dfs$p),]
                dfs$sample <- as.numeric(dfs$sample)
                dfs$F <- as.numeric(dfs$F)
                dfs$R2 <- as.numeric(dfs$R2)
                dfs$adj <- p.adjust(dfs$p, "BH")
                output$stats <- DT::renderDataTable(dfs, caption=input$cv)
                values$stats <- dfs
                
                ## Compact tree section
                ## Show (compact) trees if significant
                sigdfs <- dfs[dfs$adj < 0.05,]
                if (dim(sigdfs)[1]!=0) {
                    sigIDs <- row.names(sigdfs)
                    output$showSigTree <- renderPlot({
                        wrap_plots(lapply(sigIDs, function(id) {
                            tmpmet <- metals[[id]]
                            tmpmet <- cbind(row.names(tmpmet), tmpmet)
                            ggtree(trels[[id]], layout=input$layout,
                                   branch.length="none") %<+% tmpmet +
                                geom_tippoint(aes_string(color=input$cv), size=input$size)+
                                ggtitle(id)+theme(legend.position="bottom")+
                                scale_color_viridis_d(option="F")
                        }), nrow=1) + 
                        plot_layout(guides = "collect") & theme(legend.position = 'bottom')
                    }, ## End of plotting panel
                    bg = 'white',
                    res=96,
                    height=input$height,
                    width=input$width)
                }
                updateTabsetPanel(session, "main_tab", selected = "Inspect")
            } else {
                ## Make heatmap based on F-values (or other statistics) for multiple CVs
                if (input$hmcell=="pseudoF") {
                    FMat <- do.call(rbind, lapply(ps, function(x) x[["F"]]))
                    HMNAME <- "Fvalue"
                } else if (input$hmcell=="P") {
                    FMat <- do.call(rbind, lapply(ps, function(x) x[["p"]]))
                    HMNAME <- "P"
                } else {
                    FMat <- do.call(rbind, lapply(ps, function(x) x[["r2"]]))
                    HMNAME <- "R2"
                }
                FMat <- FMat |> `row.names<-`(lapply(ps,
                　　function(x) paste0(x[["sp"]], " (",x[["totdf"]],")") |>
                    unlist()))
                colnames(FMat) <- input$cv

                ## Remove all NA row and col
                FMat <- FMat[rowSums(is.na(FMat))!=ncol(FMat),]
                FMat <- FMat |> as.data.frame() |>
                    dplyr::select(which(colSums(is.na(FMat))!=nrow(FMat)))

                ## [SAVING part, omitted]
                # prefixdate <- format(Sys.time(), "%b%d%X%Y") %>%  gsub(":","",.)
                # save(file=paste0("FMat",prefixdate,".rda"), FMat)
                values$mat <- FMat
                FHeatmap <- Heatmap(FMat,
                                    name=HMNAME, border=FALSE,
                                    rect_gp = gpar(col = "grey20", lwd = 2),
                                    column_names_side="top",
                                    row_names_max_width = unit(10, "cm"))
                values$heatmap <- FHeatmap
                tryCatch(
                    {
                        output$showGenePlot <- renderPlot({FHeatmap},
                                                          bg = 'white', res=96, height=input$height,
                                                          width=ifelse(input$width>800, input$width, 800))
                    },
                    error=function (e) {print(e)}
                )

                ## Move to CNV tab
                updateTabsetPanel(session, "main_tab", selected = "CNVs")
            }
        } else {
            updateLog("No comparison made")
        }

        waiter_hide()
    }) ## End of alladonis

    ## Select species and input IDs in sp_id
    observeEvent(input$species, {
        values$sp_id <- values$rev_names[input$species]
        values$species <- input$species
        updateLog(paste0("Selecting ", input$species))
    })

    ## Tree showing
    observeEvent(input$show,{
        waiter_show(html = spin_3(), color = "white")
        reset(output)
        if( is.null(input$cv) ) {waiter_hide();return(1)}

        if (input$use_geom_fruit) {
            cat_subtle("# Using ggtreeExtra::geom_fruit\n")
            show_cv <- input$cv
        } else {
            show_cv <- input$cv[1]
        }
        shinyjs::runjs("window.scrollTo(0, 0)")

        sp_id <- values$sp_id

        cat_subtle("# Showing ", sp_id, " ", show_cv, "\n", sep="")

        if (input$all>1) {## Whole species
            ## [NOTE] Adonis will not be performed.
            ## We cannot do multi-meta because of the small plot produced.
            if (input$use_geom_fruit) {return(1);waiter_hide();}

            plot_list <- list()
            for (i in values$ordered[seq_len(input$all)]) {
                tre <- stana@treeList[[i]]
                if (is.null(tre$tip.label) | length(tre$tip.label)==1) {
                    ## DO NOTHING
                } else {
                    meta <- values$meta

                    for (kk in tre$tip.label) {
                        if (!kk %in% row.names(meta)) {
                            cur <- row.names(meta)
                            meta <- rbind(meta, rep(NA, ncol(meta)))
                            row.names(meta) <- c(cur, kk)
                        }
                    }

                    meta <- meta[tre$tip.label,]

                    ## Layout for `all` type should be fixed?
                    # meta$label <- NULL
                    # meta <- cbind(meta[, 1], meta)
                    p <- ggtree(tre, layout=input$layout,
                                branch.length=ifelse(input$cladogram,
                                                     "none", "branch.length")) %<+% meta +
                        geom_tippoint(aes_string(color=show_cv), size=input$size)+
                        ggtitle(i)+theme(legend.position="bottom") ## Title wil be just ID in MIDASDB
                    if (input[[show_cv]]) {
                        p <- p + scale_color_gradient(low="blue",
                                                      high="red")
                    } else {
                        p <- p + scale_color_brewer(palette=scicoPal)
                    }
                    plot_list[[i]] <- p
                }
            }
            panelShow <- patchwork::wrap_plots(plot_list, ncol=input$colnum)
            output$showPlot <- renderPlot({panelShow}, bg = 'white', res=96, height=input$height, width=input$width)
            waiter_hide()
        } else {## One species
            # https://mastering-shiny.org/action-dynamic.html
            ## Show title of species
            updateLog(paste0("Showing tree for ", input$species))
            
            values$species <- input$species
            output$species_name <- renderText(values$species)
            values$sp_id <- sp_id
            ## Load the tree
            ## If tree is not available, use KO (or gene)-NMF approach
            if (is.null(values$stana@treeList[[sp_id]]) | input$use_nmf) {
                ko_df_filt <- values$stana@kos[[sp_id]]
                if (dim(ko_df_filt)[1]<=1) {
                    shiny::showNotification("Too small sample");
                    waiter_hide();
                    return(1)                	
                }
                if (!nnlm) {
				    ko_df_filt <- ko_df_filt[apply(ko_df_filt, 1, function(x) sum(x)!=0),]
				    ko_df_filt <- ko_df_filt[,apply(ko_df_filt, 2, function(x) sum(x)!=0)]                	
                }

                output$showPlot <- NULL
                output$adonis <- NULL

                cat_subtle("# Using gene or KO table to calculate factor matrix\n")
                cat_subtle("# Will take long time depending on the input matrix dimension\n")
                if (nnlm) {
                	if (input$rank==0) {
    	                estimate_range <- 1:6
		                cat_subtle("# estimate_range: ", paste0(estimate_range, collapse="/"), "\n", sep="")
		                ## NA freq is fixed to 0.3
		                ## Just one run
			    		A <- as.matrix(ko_df_filt)
						already_na <- which(is.na(A))
						allind <- seq_len(length(A))
						newind <- allind[!(allind %in% already_na)]
						ind <- sample(newind, 0.3 * length(newind));
						A2 <- A;
						A2[ind] <- NA;

						err <- sapply(X = estimate_range,
						              FUN = function(k) {
						                z <- nnmf(A2, k);
						                c(mean((with(z, W %*% H)[ind] - A[ind])^2), tail(z$mse, 1));
						              }
						);
						rank <- which.min(err[1,])	
				        cat_subtle("# Chosen rank: ", rank, "\n", sep="")
                	} else {
                		rank <- input$rank
                	}
                	z <- nnmf(as.matrix(ko_df_filt), rank);
                	coefMat <- z$H
		            relab <- apply(coefMat, 2, function(x) x / sum(x))	

                } else {
                	if (input$rank==0) {
    	                estimate_range <- 1:6
		                cat_subtle("# estimate_range: ", paste0(estimate_range, collapse="/"), "\n", sep="")
			            tryCatch(
			                {
		                    	test <- nmfEstimateRank(as.matrix(ko_df_filt),
		                    		range=estimate_range, method="snmf/r")   	
			                },
			                error=function (e) {print(e)}
			            )
			            
			            if (!exists("test")) {
			            	## Loose checking
			            	waiter_hide();
			            	return(1);
			            }

		            	val <- test$measures[, "cophenetic"]
		                b <- -1
				        for (i in seq_along(val)) {
				            if (is.na(val[i])) {
				                next
				            } else {
				                if (val[i] > b) {
				                    b <- val[i]
				                } else {
				                break
				                }
				            }
				        }
		         	    rank <- estimate_range[i]
				        cat_subtle("# Chosen rank: ", rank, "\n", sep="")
                	} else {
                		rank <- input$rank
                	}
		            res <- NMF::nmf(ko_df_filt, rank = rank, seed = 1, method="snmf/r")
		            coefMat <- coef(res)
		            relab <- apply(coef(res), 2, function(x) x / sum(x))	
                }
                
	            ## Metadata processing
                meta <- values$meta

	            ## Add NA row to unknown label in tree
	            ## not in metadata
	            for (i in colnames(relab)) {
	                if (!i %in% row.names(meta)) {
	                    cur <- row.names(meta)
	                    meta <- rbind(meta, rep(NA, ncol(meta)))
	                    row.names(meta) <- c(cur, i)
	                }
	            }

	            for (tmp_show_cv in show_cv) {
	                change_cv <- meta[[tmp_show_cv]]
	                change_cv[change_cv==""] <- NA
	                change_cv[change_cv=="-"] <- NA
	                change_cv[is.infinite(change_cv)] <- NA

	                if (input[[tmp_show_cv]]) {
	                    meta[[tmp_show_cv]] <- as.numeric(change_cv)
	                } else {
	                    meta[[tmp_show_cv]] <- as.factor(change_cv)
	                }
	                if (is.na(meta[[tmp_show_cv]]) |> sum() == dim(meta)[1]) {
	                    cat_subtle("# One of the covariates are all NA\n")
	                    updateLog("One of the covariates are all NA, exiting");
	                    waiter_hide()
	                    return(1)
	                }
	                if (!tmp_show_cv %in% colnames(meta)) {
	                    shiny::showNotification("Selected column is not present in metadata")
	                    waiter_hide();
	                    return(1)
	                }
	            }
	            
	            
            	stb <- data.frame(t(relab))
				colnames(stb) <- as.character(seq_len(ncol(stb)))
			    stb[["sample"]] <- row.names(stb)
                cd <- meta[[show_cv]] %>% setNames(row.names(meta))
                stb[["group"]] <- cd[stb[["sample"]]]

        	    ## Not showing sample label, instead facet by group
                if (is.numeric(meta[[show_cv]])) {
                    melted <- reshape2::melt(stb, measure.vars=as.character(seq_len(rank)))
                    panelShow <- ggplot(melted, aes(x=group, y=value)) + 
                        geom_point()+
                        facet_grid(. ~ variable, scale="free")+
                        cowplot::theme_cowplot()+ cowplot::panel_border()
                } else {
                    melted <- reshape2::melt(stb)
                    panelShow <- ggplot(melted, aes(fill=variable, y=value, x=sample)) + 
                        geom_col(position="fill")+
                        facet_grid(. ~ group, scale="free")+
                        cowplot::theme_cowplot()+ cowplot::panel_border()+
                        theme(axis.text.x = element_blank())                    
                }


                values$p <- panelShow
                output$showPlot <- renderPlot({panelShow},
                              bg = 'white',
                              res=96,
                              height=input$height,
                              width=input$width)
                waiter_hide();
                return(1)
            } ## KO-NMF approach
            
            if (input$use_ko_dist) {
                cat_subtle("# Using KO table to calculate distance between samples\n")
                ## Load the precalculated KO table
                ko_df_filt <- values$stana@kos[[sp_id]]
                if (dim(ko_df_filt)[1]<=1) {
                    shiny::showNotification("Too small sample");
                    waiter_hide();
                    return(1)                   
                }
                dist_mat <- dist(t(ko_df_filt))
                tre <- phangorn::upgma(dist_mat)
            } else {
                cat_subtle("# Reading tree for ", sp_id, "\n")
                tre <- values$stana@treeList[[sp_id]]
            }

            ## Processing if tree has 0 or 1 sample
            if (is.null(tre$tip.label) | length(tre$tip.label)==1) {
                waiter_hide();
                shiny::showNotification("Not enough samples are profiled")
                nsample <- tre$tip.label |> length()
                output$showPlot <- NULL
                output$adonis <- NULL
                output$nSample <- renderText(paste0("Profiled in : ", nsample, " samples"))
                return(1);
            }

            ## Metadata processing
            meta <- values$meta

            ## Add NA row to unknown label in tree
            ## not in metadata
            for (i in tre$tip.label) {
                if (!i %in% row.names(meta)) {
                    cur <- row.names(meta)
                    meta <- rbind(meta, rep(NA, ncol(meta)))
                    row.names(meta) <- c(cur, i)
                }
            }

            for (tmp_show_cv in show_cv) {
                change_cv <- meta[[tmp_show_cv]]
                change_cv[change_cv==""] <- NA
                change_cv[change_cv=="-"] <- NA
                change_cv[is.infinite(change_cv)] <- NA

                if (input[[tmp_show_cv]]) {
                    meta[[tmp_show_cv]] <- as.numeric(change_cv)
                } else {
                    meta[[tmp_show_cv]] <- as.factor(change_cv)
                }
                if (is.na(meta[[tmp_show_cv]]) |> sum() == dim(meta)[1]) {
                    cat_subtle("# One of the covariates are all NA\n")
                    updateLog("One of the covariates are all NA, exiting");
                    waiter_hide()
                    return(1)
                }
                if (!tmp_show_cv %in% colnames(meta)) {
                    shiny::showNotification("Selected column is not present in metadata")
                    waiter_hide();
                    return(1)
                }
            }

            meta <- meta[tre$tip.label,]

            ## Test PERMANOVA
            if (!input$use_ko_dist) {
                if (input$use_pairwise) {
                    cat_subtle("# Pairwise mode enabled\n")
                    wholeLevel <- unique(meta[[show_cv[[1]]]])
                    wholeLevel <- wholeLevel[!is.na(wholeLevel)]
                    cat_subtle("# ", length(wholeLevel), " level present\n", sep="")
                    comps <- combn(wholeLevel, m=2)
                    aresall <- do.call(rbind, lapply(seq_len(ncol(comps)), function(l) {
                        levs <- comps[,l]
                        levs1 <- levs[1]
                        levs2 <- levs[2]
                        
                        inc_samples <- row.names(subset(meta, meta[[show_cv[[1]]]] %in% levs))
                        ## Multiple subtree?
                        # subt <- castor::get_subtree_with_tips(tre, inc_samples)$subtree
                        dist_mat <- as.dist(ape::cophenetic.phylo(tre))
                        dist_mat <- as.dist(as.matrix(dist_mat)[inc_samples, inc_samples])
                        subset_meta2 <- subset(meta, meta[[show_cv[[1]]]] %in% levs)
                        subset_meta2 <- subset_meta2[inc_samples, ]
                        
                        if ((subset_meta2[[show_cv[1]]] |> na.omit() |> unique() |> length()) != 1) {
                            adonis_res <- vegan::adonis2(as.formula(paste0("dist_mat ~ ", show_cv[[1]])),
                                                         data=subset_meta2, na.action = na.omit)
                        } else {
                            return(NULL)
                        }
                        ado_df <- data.frame(adonis_res, check.names=F)
                        ado_df[["comparison"]] <- paste0(show_cv[[1]],":",levs1,"_vs_",levs2)
                        return(ado_df)
                    }))
                    cat_subtle("# adonis finished\n")
                    if (!is.null(aresall)) {
                        ado_df <- data.frame(aresall, check.names=F)
                        output$adonis <- DT::renderDataTable(ado_df)
                    } else {
                        cat_subtle("# No adonis results produced\n")
                    }
                } else { ## End of pairwise
                    dist_mat <- as.dist(ape::cophenetic.phylo(tre))
                }
            }
            
            if (!input$use_pairwise) {
                subset_meta <- meta[tre$tip.label,]
                
                ## [CAUTION] If NA is present in metadata, they are omitted!
                ## [CAUTION] In the `show` mode, test will be performed only on the first covariate
                ## [CAUTION] If Pair-wise, first the tree will be subset to the ID, and the tests
                ## were performed for each pair.
                if ((subset_meta[[show_cv[1]]] |> na.omit() |> unique() |> length()) != 1) {
                    adonis_res <- vegan::adonis2(as.formula(paste0("dist_mat ~ ", show_cv[[1]])),
                                                 data=subset_meta, na.action = na.omit)
                } else {
                    adonis_res <- NULL
                }
                cat_subtle("# adonis finished\n")
                if (!is.null(adonis_res)) {
                    ado_df <- data.frame(adonis_res, check.names=F)
                    output$adonis <- DT::renderDataTable(ado_df)
                }
            }

            if (input$cladogram) {
                branch.length = "none"
            } else {
                branch.length = "branch.length"
            }
            if (input$use_geom_fruit) {
                ## Can plot multiple geoms
                if (input$layout=="rectangular") {
                    cat_subtle("# Force setting layout to 'circular'\n")
                    flyt <- "circular"
                } else {
                    flyt <- input$layout
                }

                ## Select random palette from scico
                scp <- sample(RColorBrewer::brewer.pal.info %>% row.names(),
                    length(show_cv), replace=FALSE)
                names(scp) <- show_cv

                ## Select random shape from ggstar
                starsh <- sample(1:30, length(show_cv), replace=FALSE)
                names(starsh) <- show_cv

                p <- ggtree(tre, layout=flyt, branch.length=branch.length)
                subset_meta[["id"]] <- row.names(subset_meta)
                for (tmp_show_cv in show_cv){
                    cat_subtle("# Processing ", tmp_show_cv, "\n", sep="")
                    if (is.numeric(subset_meta[[tmp_show_cv]])) {
                        cat_subtle("# Numeric\n")
                        ## If is numeric
                        p <- p + geom_fruit(
                            data = subset_meta,
                            geom = geom_col,
                            mapping = aes(y=id, fill=.data[[tmp_show_cv]], x=.data[[tmp_show_cv]]),
                        ) +
                            scale_fill_distiller(palette = scp[tmp_show_cv]) +
                            new_scale_fill()
                    } else {
                        cat_subtle("# Discrete\n")
                        ## If is discrete
                        p <- p + geom_fruit(
                            data = subset_meta,
                            geom = geom_star,
                            size = input$size,
                            mapping = aes(y=id, fill=.data[[tmp_show_cv]]),
                            starshape = starsh[tmp_show_cv]
                        ) +
                            scale_fill_brewer(palette = scp[tmp_show_cv]) +
                            new_scale_fill()
                    }
                }
            } else {
                # meta$label <- NULL
                # meta <- cbind(meta[, 1], meta)
                p <- ggtree(tre, layout=input$layout, branch.length=branch.length) %<+% meta +
                    geom_tippoint(aes_string(color=show_cv[1]), size=input$size)
                if (is.numeric(meta[[show_cv[1]]])) {
                    p <- p + scale_color_gradient(low="blue",
                                                  high="red")
                } else {
                    p <- p + scale_color_brewer(palette = pal)
                }
            }

            showCol <- colnames(meta)
            showCol <- showCol[showCol!="label"]
            ## Append species name as title
            p <- p + ggtitle(values$species)

            values$p <- p#p@net+ theme(plot.background = element_rect(fill = "transparent",colour = NA))
            panelShow <- p#p@net+ theme(plot.background = element_rect(fill = "transparent",colour = NA))
            output$showPlot <- renderPlot({panelShow},
                                          bg = 'white',
                                          res=96,
                                          height=input$height,
                                          width=input$width)
            nsample <- tre$tip.label |> length()
            output$nSample <- renderText(paste0("Profiled in : ", nsample, " samples"))
            waiter_hide()
        }
    }
    )

    ## Gene part (DEA, Heatmap)
    observeEvent(input$genes,{
        waiter_show(html = spin_3(), color = "white")
        ## Load the precalculated KO table
        show_cv <- input$cv[1]
        if (length(show_cv)==0) {waiter_hide();return(1)}
        ## No KO file available
        updateTabsetPanel(session, "main_tab", selected = "CNVs")
        ko_df_filt <- values$stana@kos[[values$sp_id]]
        if(!startsWith(row.names(ko_df_filt)[1], "ko")) {showNotification("Seems like not KO table");
            waiter_hide();return(1);}
        if (dim(ko_df_filt)[1]==0) {waiter_hide();return(1)}
        ## Metadata processing in GENE

        meta <- values$meta
        whole_meta_num <- nrow(meta)
        change_cv <- meta[[show_cv]]
        change_cv[change_cv==""] <- NA
        change_cv[change_cv=="-"] <- NA
        change_cv[is.infinite(change_cv)] <- NA

        if (input[[show_cv]]) {
            meta[[show_cv]] <- as.numeric(change_cv)
        } else {
            meta[[show_cv]] <- as.factor(change_cv)
        }

        meta <- meta[intersect(row.names(meta), colnames(ko_df_filt)),]


        if (is.na(meta[[show_cv]]) |> sum() == dim(meta)[1]) {
            cat_subtle("Covariates are all NA\n")
            updateLog("Covariates are all NA, exiting");
            waiter_hide()
            return(1)
        }

        if (!show_cv %in% colnames(meta)) {
            shiny::showNotification("Selected column is not present in metadata")
            waiter_hide();
            return(1)
        }

        ## Subset to metadata-included samples
        gene_meta <- meta[intersect(row.names(meta),
                                    colnames(ko_df_filt)),]
        gene_meta <- gene_meta[!is.na(gene_meta[[show_cv]]),]
        sel_cv <- gene_meta[[show_cv]]
        cur_cv <- show_cv
        showNum <- ncol(ko_df_filt)

        ## If factor, change based on names
        if (is.factor(sel_cv)) {
            if (input$lfc & (length(unique(sel_cv))!=2)) {
                showNotification("T-statistic is available for the group with two levels");waiter_hide();return(1);
            }
            change_dic <- gene_meta[[show_cv]] |> setNames(gene_meta[,1])
            spl <- change_dic[colnames(ko_df_filt)] |> na.omit()
        } else {
            ## If numeric and two levels...
            if (length(unique(sel_cv))==2) {
                change_dic <- gene_meta[[show_cv]] |> setNames(gene_meta[,1])
                spl <- change_dic[colnames(ko_df_filt)] |> na.omit()
            } else {
                ## We divide numeric by median
                cat_subtle("# Dividing by median ... ")
                med_val <- median(gene_meta[[show_cv]], na.rm=TRUE) ## although already removed

                cat_subtle("# ", med_val, " in ", dim(gene_meta)[1],"\n", sep="")
                change_dic <- gene_meta[[show_cv]] > med_val
                change_dic <- ifelse(change_dic, "High","Low") |>
                    setNames(gene_meta[,1])
                sel_cv <- change_dic
                cur_cv <- paste0(show_cv, "_binned")
                gene_meta[[cur_cv]] <- change_dic
                spl <- change_dic[colnames(ko_df_filt)] |> na.omit()

            }
        }
        ## Prepare for the BG set, should be modified if specified in certain environments
        kopgsea <- kop
        kopgsea <- kopgsea |> filter(startsWith(V1, "path:ko"))
        kopgsea$V1 <- kopgsea$V1 |> strsplit(":") |>
            vapply("[",2,FUN.VALUE="a")
        
        paths <- unique(kopgsea$V1)
        nl <- lapply(paths, function(p) {
            tmp <- subset(kopgsea, kopgsea$V1 == p)
            tmp$V2
        })
        names(nl) <- paths
        if (input$lfc) {
            ## Calculate LFC if two variables!
            a <- unique(sel_cv)[1]
            b <- unique(sel_cv)[2]
            comp_label <- paste0(show_cv, ": ", a, " / ", b)

            aa <- subset(gene_meta, gene_meta[[cur_cv]]==a) |> row.names()
            bb <- subset(gene_meta, gene_meta[[cur_cv]]==b) |> row.names()

            if (length(aa)==0) {waiter_hide();return(1)}
            if (length(bb)==0) {waiter_hide();return(1)}

            if (!input$gsea) {
                cat_subtle("# Calculating Wilcoxon rank-sum tests\n")
                pvalAdjustMethod <- ifelse(input$pvalAdjustMethod, "BH", "none")
                apvs <- p.adjust(apply(ko_df_filt, 1, function(x) {
                    wt <- exactRankTests::wilcox.exact(x[aa], x[bb]); wt$p.value}),
                    pvalAdjustMethod)
                cat_subtle("# Significant KOs: ", length(apvs[apvs<0.05]), "\n", sep="")
                if (length(apvs[apvs<0.05])==0) {
                    waiter_hide()
                    reset_enr(output)
                    return(1);showNotification("No significant KOs");
                }
            }
            ## eps values
            ko_sum <- mod.T(ko_df_filt, aa, bb)
            if (is.null(ko_sum)) {waiter_hide(); return(1); showNotification("Matrix is empty.")}
            ko_sum <- ko_sum[order(ko_sum, decreasing=TRUE)]
            

            if (input$gsea) {
                updateLog(paste0("Performing GSEA for ", values$sp_id))




                res <- fgsea::fgsea(nl, ko_sum, nproc=1)
                enr <- data.frame(res)

                ## Changing list to char
                enr[[8]] <- unlist(lapply(enr[[8]], function(x) paste0(x, collapse="/")))
                enr <- enr[order(enr$padj), ]

            } else {
                updateLog(paste0("Performing ORA for ", values$sp_id))


                if (input$updown) {
                    candint <- intersect(names(ko_sum[ko_sum>0]), names(apvs[apvs < 0.05]))
                } else {
                    candint <- intersect(names(ko_sum[ko_sum<0]), names(apvs[apvs < 0.05]))
                }
                cat_subtle("Length of candint: ", length(candint), "\n", sep="")
                if (length(candint)==0) {waiter_hide();return(1);}
                
                enr <- fgsea::fora(nl, candint, universe=names(ko_sum))

                res <- data.frame(enr)
                if (is.null(enr)) {waiter_hide();return(1);}
            }
        } else {
            comp_label <- paste0("Sum of all samples")
            ## Sum values for all the variables
            ko_sum <- apply(ko_df_filt, 1, sum)
            top50 <- names(ko_sum[order(ko_sum, decreasing=TRUE)])
            ## Performing ORA on all KOs
            cat_subtle("# Performing EA (for all the KOs)\n")
            enr <- fgsea::fora(nl, top50, universe=names(ko_sum))
            res <- data.frame(enr)
        }

        cat_subtle("# Inserting named vector of values to reactive\n")

        pnames <- nmls$V2 |> setNames(nmls$V1)
        enr$Description <- pnames[enr$pathway]
        enr <- enr %>% relocate(Description)

        values$lfc <- ko_sum
        values$enr <- enr
        output$showEnr <- DT::renderDataTable(enr, caption=comp_label)
        all_path <- enr$Description
        updateSelectInput(session, "path_selector",
                          choices = all_path)

        ## Visualize heatmap
        ## K-means for wordclouds

        ## Reduce plot (variable features)
        if (input$reduce!=0) {
          if (input$reduce<100) {red <- 100} else {red <- input$reduce}
          top100 <- matrixStats::rowVars(ko_df_filt) %>% 
            sort(decreasing=TRUE) %>% 
            head(red) %>% names()
          
          ko_df_filt <- ko_df_filt[top100, ]
        }
        
        km = kmeans(ko_df_filt, centers = input$kmeans)$cluster
        gene_list <- split(row.names(ko_df_filt), km)
        gene_list <- lapply(gene_list, function(x) {
            x[!is.na(x)]
        })

        annotList <- list()
        for (i in names(gene_list)) {
            maps <- (kop |> dplyr::filter(V2 %in% gene_list[[i]]))$V1
            all_maps <- maps[grepl("ko", maps)] |> sort() |> unique()
            paths <- strsplit(all_maps, "path:") |> sapply("[", 2)
            annotList[[i]] <-   (nmls |> dplyr::filter(V1 %in% paths))$V2
        }
        hm <- Heatmap(ko_df_filt[, names(spl)],
            row_split = km,
            column_split=spl,
            use_raster=TRUE,
            raster_by_magick=TRUE,
            show_column_names = FALSE,
            border=TRUE,
            name="Copy number",
            show_row_names = FALSE,
            show_row_dend = FALSE) +
            rowAnnotation(
                keywords = simplifyEnrichment::anno_word_cloud(align_to = km,
                    term=annotList,
                    exclude_words=c(
                        "pathway",
                        "metabolism",
                        "signaling",
                        "biosynthesis",
                        "degradation"),
                    fontsize_range=c(10,20),
                    max_words = 30)
            )
        output$showGenePlot <- renderPlot({draw(hm, column_title=show_cv)},
            bg = 'white',
            res=72,
            height=input$height,
            width=ifelse(input$width>800,input$width,800))
        values$heatmap <- hm
        waiter_hide()})

    ## All GSEA (output mat)
    observeEvent(input$allgsea, {
        if (!input$gsea) {showNotification("Please turn on GSEA option");return(1)}
        if (!input$lfc) {showNotification("Please turn on T-stat option");return(1)}

        waiter_show(html = spin_3(), color = "white")
        ## Load the precalculated KO table
        list_of_gsea <- list()
        level1 <- NULL
        level2 <- NULL
        for (sp in values$current_list_of_species) {
            sp_id <- values$rev_names[sp]
            ## Load the tree
            tre <- values$stana@treeList[[sp_id]]
            labnum <- length(tre$tip.label)

            ## Metadata processing in GENE

            meta <- values$meta
            whole_meta_num <- nrow(meta)

            if (input$thresh>=1) {
                curthre <- input$thresh
            } else {
                curthre <- whole_meta_num * input$thresh
            }

            if (labnum > curthre) {

                show_cv <- input$cv[1]
                ## No KO file available
                ko_df_filt <- values$stana@kos[[sp_id]]



                meta <- values$meta

                change_cv <- meta[[show_cv]]
                change_cv[change_cv==""] <- NA
                change_cv[change_cv=="-"] <- NA
                change_cv[is.infinite(change_cv)] <- NA


                if (input[[show_cv]]) {
                    meta[[show_cv]] <- as.numeric(change_cv)
                } else {
                    meta[[show_cv]] <- as.factor(change_cv)
                }

                if (is.na(meta[[show_cv]]) |> sum() == dim(meta)[1]) {
                    cat_subtle("Covariates are all NA\n")
                    updateLog("Covariates are all NA, exiting");
                    # waiter_hide()
                    next
                }

                if (!show_cv %in% colnames(meta)) {
                    shiny::showNotification("Selected column is not present in metadata")
                    next
                }

                ## Subset to metadata-included samples
                gene_meta <- meta[intersect(row.names(meta),
                                            colnames(ko_df_filt)),]
                gene_meta <- gene_meta[!is.na(gene_meta[[show_cv]]),]
                sel_cv <- gene_meta[[show_cv]]
                cur_cv <- show_cv
                showNum <- ncol(ko_df_filt)

                ## If factor, change based on names
                if (is.factor(sel_cv)) {
                    if (input$lfc & (length(unique(sel_cv))!=2)) {
                        showNotification("Comparison in GSEA is available for the group with two levels");
                        waiter_hide();
                        break;
                    }
                    change_dic <- gene_meta[[show_cv]] |> setNames(gene_meta[,1])
                    spl <- change_dic[colnames(ko_df_filt)] |> na.omit()
                } else {
                    ## If numeric and two levels...
                    if (length(unique(sel_cv))==2) {
                        change_dic <- gene_meta[[show_cv]] |> setNames(gene_meta[,1])
                        spl <- change_dic[colnames(ko_df_filt)] |> na.omit()
                    } else {
                        ## We divide numeric by median
                        cat_subtle("# Dividing by median ... ")
                        med_val <- median(gene_meta[[show_cv]], na.rm=TRUE) ## although already removed

                        cat_subtle("# ", med_val, " in ", dim(gene_meta)[1], "\n", sep="")
                        change_dic <- gene_meta[[show_cv]] > med_val
                        change_dic <- ifelse(change_dic, "High","Low") |>
                            setNames(gene_meta[,1])
                        sel_cv <- change_dic
                        cur_cv <- paste0(show_cv, "_binned")
                        gene_meta[[cur_cv]] <- change_dic
                        spl <- change_dic[colnames(ko_df_filt)] |> na.omit()

                    }
                }
                ## Prepare for the BG set, should be modified if specified in certain environments
                kopgsea <- kop
                kopgsea <- kopgsea |> filter(startsWith(V1, "path:ko"))
                kopgsea$V1 <- kopgsea$V1 |> strsplit(":") |>
                    vapply("[",2,FUN.VALUE="a")
                
                paths <- unique(kopgsea$V1)
                nl <- lapply(paths, function(p) {
                    tmp <- subset(kopgsea, kopgsea$V1 == p)
                    tmp$V2
                })
                names(nl) <- paths
                if (input$lfc) {
                    ## Calculate LFC if two variables!
                    ## Fix levels for all the comparison
                    if (is.null(level1)) {
                        a <- unique(sel_cv)[1]
                        level1 <- TRUE
                    }
                    if (is.null(level2)) {
                        b <- unique(sel_cv)[2]
                        level2 <- TRUE
                    }

                    comp_label <- paste0(show_cv, ": ", a, " / ", b)
                    cat_subtle("# ", comp_label, "\n", sep="")
                    aa <- subset(gene_meta, gene_meta[[cur_cv]]==a) |> row.names()
                    bb <- subset(gene_meta, gene_meta[[cur_cv]]==b) |> row.names()

                    if (length(aa)==0) {next;}
                    if (length(bb)==0) {next;}
                    comp_label <- paste0(show_cv, ": ", a, " / ", b)
                    
                    ko_sum <- mod.T(ko_df_filt, aa, bb)
                    if (is.null(ko_sum)) {next; showNotification("Matrix is empty.")}

                    ## Perform GSEA (it will take time)?
                    ko_sum <- ko_sum[order(ko_sum, decreasing=TRUE)]
                    # top50 <- names(abs(ko_sum)[order(abs(ko_sum), decreasing=TRUE)]) |> strsplit(":") |> sapply("[",2)
                    updateLog(paste0("Performing GSEA for ", sp_id))
                    cat_subtle("# GSEAMat: Cleared the threshold (", curthre, "), performing GSEA for ", sp_id, " (", labnum, "\n", sep="")
                

                    if (length(ko_sum)==0) {next;}
                    withProgress(
                        res <- fgsea::fgsea(nl, ko_sum, nproc=1),
                        message=paste0("Performing GSEA: ", sp_id)
                    )
                    
                    enr <- data.frame(res)

                    ## Changing list to char
                    enr[[8]] <- unlist(lapply(enr[[8]], function(x) paste0(x, collapse="/")))
                    ## Saving function
                    # save(file=paste0(input$site,"_",cur_cv,"_",sp,".rda"), enr)
                } else {
                    comp_label <- paste0("Sum of all samples")
                    ## Sum values for all the variables
                    ko_sum <- apply(ko_df_filt, 1, sum)
                    top50 <- names(ko_sum[order(ko_sum, decreasing=TRUE)]) |>
                        strsplit(":") |> sapply("[", 2)
                    ## Performing ORA on all KOs
                    cat_subtle("# Performing EA\n")
                    enr <- fgsea::fora(nl, top50, names(ko_sum)) %>% data.frame()
                }
                list_of_gsea[[sp]] <- enr
            } else {
                cat_subtle("# GSEAMat: Not cleared the threshold (", curthre, "), ", sp_id," (", labnum, "\n", sep="")
            }
        }
        ## allgsea plot
        if (length(list_of_gsea)!=0) {
            updateTabsetPanel(session, "main_tab", selected = "CNVs")
            tryCatch(
                {
                    GSEA_ALL_HEATMAP <- return_all_gsea_heatmap(list_of_gsea, sel="NES", title=comp_label)
                    output$showGenePlot <- renderPlot({GSEA_ALL_HEATMAP},
                                                      bg = 'white', res=96,
                                                      height=ifelse(input$height>2000, input$height, 2000),
                                                      width=ifelse(input$width>1500, input$width, 1500))
                    values$heatmap <- GSEA_ALL_HEATMAP
                    waiter_hide()
                    values$gseamat <- return_all_gsea_heatmap(list_of_gsea, sel="NES",retTable=TRUE)
                    output$downloadGSEAMat <- downloadHandler(
                        filename = function(){"gsea_mat.csv"},
                        content = function(fname){
                            write.csv(values$gseamat, fname)
                        }
                    )
                },
                error=function (e) {print(e)}
            )
            ## Saving function
            # save(file=paste0("gsea_all_",input$site,".rda"), list_of_gsea)
        }
        waiter_hide()
    })

    ## Plotting KEGG
    observeEvent(input$kegg,{
        cat_subtle("# KEGG mode started ...\n")
        if (is.null(values$enr)) {return(1)}
        waiter_show(html = spin_3(), color = "white")
        updateTabsetPanel(session, "main_tab",
                          selected = "KEGG PATHWAY")
        cat_subtle("# ", input$path_selector, "\n", sep="")
        enr_res <- values$enr
        if (startsWith(enr_res$pathway[1], "ko")) {
            enr_res$pathway <- gsub("ko","map",enr_res$pathway)
        }
        tmp <- enr_res[ enr_res$Description==input$path_selector, ]

        if (startsWith(tmp$pathway, "map012") | startsWith(tmp$pathway, "map011")) {## If global map
            if (input$ggraph) {
                path <- pathway(gsub("map","ko",tmp$pathway)) |> process_reaction()
            } else {
                path <- pathway(gsub("map","ko",tmp$pathway)) |> process_line()
            }
            path <- path |> activate(edges) |> mutate(num = edge_numeric(values$lfc))
            dd <- (path |> activate(edges) |> data.frame())

            ## OR, use ggraph to plot the global map and corresponding LFCs
            if (input$ggraph) {
                path <- path |> activate(nodes) |>
                    mutate(x=NULL, y=NULL, ko = convert_id("ko"),
                           num=node_numeric(values$lfc),
                           comp=convert_id("compound",
                                           first_arg_comma = FALSE))
                graph <- path |> ggraph(layout="nicely")+
                    geom_node_point(aes(filter=type=="compound"))+
                    geom_edge_parallel(
                        aes(color=num),
                        end_cap=circle(1,"mm"),
                        start_cap=circle(1,"mm"),
                        arrow=arrow(length=unit(1,"mm"),type="closed"))+
                    geom_node_text(aes(label=comp,filter=type=="compound"),
                                   size=1.5,
                                   repel=TRUE,
                                   bg.colour="white")+
                    scale_edge_color_gradient(low="blue",
                                              high="red",name="Statistics")+
                    theme_graph()

                stat <- path |> activate(nodes) |> data.frame()
                statdf <- stat[,c("name","ko","num")]
                statdf <- statdf[!duplicated(statdf[,c("ko","num")]),]
            } else {
                path <- path |> activate(nodes) |>
                    mutate(ko=convert_id("ko"), num=node_numeric(values$lfc))

                graph <- path |>
                    ggraph(layout="manual", x=x, y=y) +
                    overlay_raw_map()+
                    geom_edge_link0(aes(color=num))+
                    scale_edge_color_gradient(low="blue",high="red",name="Statistics")+
                    theme_void()
                stat <- path |> activate(nodes) |> data.frame()
                statdf <- stat[,c("name","ko","num")]
                statdf <- statdf[!duplicated(statdf[,c("ko","num")]),]
            }
        } else {## If not global map
            path <- pathway(gsub("map","ko",tmp$pathway))
            path <- path |> activate(nodes) |> mutate(num = node_numeric(values$lfc),ko=convert_id("ko"))
            graph <- ggraph(path, layout="manual", x=x, y=y) +
                overlay_raw_map()+
                geom_node_rect(aes(fill=num, filter=type=="ortholog"),
                               color="black")+
                geom_node_text(aes(label=strsplit(graphics_name, ",") %>%
                                       sapply("[", 1) %>%
                                       strsplit("\\.") %>%
                                       sapply("[", 1),
                                   filter=type=="ortholog"), size=2)+
                scale_fill_gradient(low="blue",high="red",name="Statistics") +
                theme_void()
            stat <- path |> activate(nodes) |> data.frame()
            statdf <- stat[,c("name","graphics_name","ko","num")]
        }
        statdf <- statdf[!is.na(statdf$num),]
        output$KOSTAT <- DT::renderDataTable(statdf)

        # if (startsWith(tmp$ID, "map011")) {print("hoge")}
        graph <- graph + ggtitle(values$species)
        values$keggp <- graph
        ## Saving function
        # save(file="graph.rda", graph)
        output$showKEGGPlot <- renderPlot({graph},
                                          bg = 'white', res=96, height=input$height,
                                          width=ifelse(input$width>800,input$width,800))
        waiter_hide()
    })

    output$downloadKEGGPathwayPlot <- downloadHandler(
        filename = function() {
            time <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
            paste0(time,".png")
        },
        content = function(file) {
            ggsave(values$keggp, dpi=300, width=input$save_width, height=input$save_height, filename = file)
        }
    )

    output$downloadPlot <- downloadHandler(
        filename = function() {
            time <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
            paste0(time,".png")
        },
        content = function(file) {
            ggsave(values$p, dpi=300, width=input$save_width,
                   height=input$save_height, filename = file)
        }
    )

    output$downloadHeatmap <- downloadHandler(
        filename = function() {
            time <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
            paste0(time,".png")
        },
        content = function(file) {
            png(file, width=input$save_width, height=input$save_height, units="in", res=300)
            values$heatmap
            dev.off()
        }
    )

    output$downloadMat <- downloadHandler(
        filename = function() {
            time <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
            paste0(time,".tsv")
        },
        content = function(file) {
            write.table(values$mat, file=file, sep="\t", quote=FALSE)
        }
    )
    
    output$downloadStatTSV <- downloadHandler(
        filename = function() {
            time <- gsub(":", "-", gsub(" ", "-", as.character(Sys.time())))
            paste0(time,".tsv")
        },
        content = function(file) {
            write.table(values$stats, file=file, sep="\t", quote=FALSE)
        }
    )

}

shinyApp(ui, server)
