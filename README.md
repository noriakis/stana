# stana
                                           
[![R-CMD-check](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/stana/actions/workflows/R-CMD-check.yaml) [![CodeFactor](https://www.codefactor.io/repository/github/noriakis/stana/badge)](https://www.codefactor.io/repository/github/noriakis/stana)

Strain-level metagenomic analysis in R. Import and analyse, plot output of the software like [MIDAS](https://github.com/snayfach/MIDAS), [MIDAS2](https://github.com/czbiohub/MIDAS2), [metaSNV v1, metaSNV v2](https://github.com/metasnv-tool/metaSNV) and [inStrain](https://github.com/MrOlm/inStrain). Will grow into the complete package.

## Installation
```r
devtools::install_github("noriakis/stana")
```

## Databases
- [midas_db_v1.2](https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md): Contains 31,007 bacterial reference genomes clustered into 5,952 species groups.
- UHGG (in MIDAS2): Contains 286,997 genomes clustered into 4,644 species (**from human stool samples**).
- GTDB (in MIDAS2): Contains 258,406 genomes clustered into 45,555 bacterial and 2,339 archaeal species.
- proGenomes2, proGenomes3

## Example analysis
The function includes filtering of species based on clinical / environmental variables interested, calling of consensus sequencing, constructing the tree, plotting the results based on groups, and functional annotations. The below example shows the analysis of the subset of `PRJEB9584`, sequenced by HiSeq 2000. Packages including [`ggtree`](https://github.com/YuLab-SMU/ggtree), [`ComplexHeatmap`](https://github.com/jokergoo/ComplexHeatmap), [`simplyfyEnrichment`](https://github.com/jokergoo/simplifyEnrichment), [`phangorn`](https://github.com/KlausVigo/phangorn), [circlize](https://github.com/jokergoo/circlize), and [`ggraph`](https://github.com/thomasp85/ggraph) are used.

<img src="https://github.com/noriakis/software/blob/main/images/stana_example.png?raw=true" width="800px">

## TODO
- Read all data in S4 class, and perform the downstream analysis. Currently the package reads the data per function.
- Add more filters and tree inferring options in `consensusSeq()`.
- Add reasonable filtering options to `loadInStrain` and `loadmetaSNV`