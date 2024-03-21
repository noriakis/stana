
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stana

<!-- badges: start -->

[![R-CMD-check](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml)
[![CodeFactor](https://www.codefactor.io/repository/github/noriakis/stana/badge)](https://www.codefactor.io/repository/github/noriakis/stana)
<!-- badges: end -->

Metagenotyping analysis in R. Import and analyse, visualize the
metagenotyping output of the software like
[MIDAS](https://github.com/snayfach/MIDAS),
[MIDAS2](https://github.com/czbiohub/MIDAS2), [metaSNV v1, metaSNV
v2](https://github.com/metasnv-tool/metaSNV) and
[inStrain](https://github.com/MrOlm/inStrain). In general the
metagenotyping software produces the allelic count information and gene
copy number tables and the package utilizes these information to analyze
the intra-species diversity.

The detailed usage is available at
<https://noriakis.github.io/software/stana>, using `bookdown`.

## Installation

Using `devtools`:

``` r
devtools::install_github("noriakis/stana")
```

## Pipeline

<img src="https://github.com/noriakis/software/blob/main/images/stana_pipeline.png?raw=true" width="800px">

## Examples

``` r
## Using example data
library(stana)
load(system.file("extdata", "sysdata.rda", package = "stana"))

stana
#> # A stana: MIDAS2
#> # Database: uhgg
#> # Loaded directory: midas2_sample_merge_uhgg
#> # Species number: 1
#> # Group info (list): Group1/Group2
#> # Loaded SNV table: 1 ID: 100003
#> # Loaded gene table: 1 ID: 100003
#> # Loaded KO table: 1 ID: 100003
#> # Size:7623472 B
#> # 
#> # SNV description
#> # A tibble: 2 × 3
#> # Groups:   group [2]
#>   group  species_id                                               n
#>   <chr>  <chr>                                                <int>
#> 1 Group1 d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacte…     4
#> 2 Group2 d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacte…     7
getID(stana)
#> [1] "100003"

## Make example metadata
samples <- getSlot(stana, "snps")[[1]] |> colnames()
metadata <- data.frame(
    row.names=samples,
    treatment=factor(sample(1:3, length(samples), replace=TRUE)),
    marker=runif(length(samples))
)


## Set metadata
stana <- setMetadata(stana, metadata)

## Call consensus sequence
## Infer and plot tree based on metadata
stana <- stana |>
  consensusSeq(argList=list(site_prev=0.95)) |>
  inferAndPlotTree(meta=c("treatment","marker"))
#> Beginning calling for 100003
#>   Site number: 5019
#>   Profiled samples: 11
#>   Included samples: 11
getFasta(stana)[[1]]
#> 11 sequences with 896 character and 625 different site patterns.
#> The states are a c g t
getTree(stana)[[1]]
#> 
#> Phylogenetic tree with 11 tips and 10 internal nodes.
#> 
#> Tip labels:
#>   ERR1711593, ERR1711594, ERR1711596, ERR1711598, ERR1711603, ERR1711605, ...
#> 
#> Rooted; includes branch lengths.
getTreePlot(stana)[[1]]
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="1800" style="display: block; margin: auto;" />

If the gene copy number table is available like in `MIDAS` series and
`inStrain`, one can compare the functional implications of these gene
contents.

## Interactive inspection

The users can inspect metagenotyping results interactively using Shiny
based on the variables such as disease conditions
(`exportInteractive()`). One can publish the results in the hosting
services for sharing the research findings.
