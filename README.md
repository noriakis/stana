
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stana

<!-- badges: start -->

[![R-CMD-check](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml)
[![CodeFactor](https://www.codefactor.io/repository/github/noriakis/stana/badge)](https://www.codefactor.io/repository/github/noriakis/stana)
<!-- badges: end -->

Metagenotyping analysis in R. Import and analyse, plot output of the
software like [MIDAS](https://github.com/snayfach/MIDAS),
[MIDAS2](https://github.com/czbiohub/MIDAS2), [metaSNV v1, metaSNV
v2](https://github.com/metasnv-tool/metaSNV) and
[inStrain](https://github.com/MrOlm/inStrain). Primarly developed for
`MIDAS` and `MIDAS2`. The documentation is available using `pkgdown` at
<https://noriakis.github.io/software/stana_pkgdown>. The detailed usage
is available [here](https://noriakis.github.io/software/stana), using
`bookdown`.

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
load("inst/extdata/sysdata.rda")
stana@ids
#> [1] "100003"
stana <- stana |>
  consensusSeq(argList=list(site_prev=0.95)) |>
  plotTree()
#> Beginning calling for 100003
#>   Site number: 5019
#>   Profiled samples: 11
#>   Included samples: 11
stana@fastaList[[1]]
#> 11 sequences with 896 character and 625 different site patterns.
#> The states are a c g t
stana@treeList[[1]]
#> 
#> Phylogenetic tree with 11 tips and 9 internal nodes.
#> 
#> Tip labels:
#>   ERR1711593, ERR1711594, ERR1711596, ERR1711598, ERR1711603, ERR1711605, ...
#> 
#> Unrooted; includes branch lengths.
```

## Interactive inspection

The users can inspect metagenotyping results interactively using Shiny
based on the variables such as disease conditions
(`exportInteractive()`).
