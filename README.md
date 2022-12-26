# stana

Strain-level metagenomic analysis in R. Import and analysize, plot output of the software like [MIDAS](https://github.com/snayfach/MIDAS), [MIDAS2](https://github.com/czbiohub/MIDAS2), metaSNV v1 and metaSNV v2. Will grow into the complete package.

## Databases
- [midas_db_v1.2](https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md): Contains 31,007 bacterial reference genomes clustered into 5,952 species groups.
- UHGG (in MIDAS2): Contains 286,997 genomes clustered into 4,644 species (**from human stool samples**).
- GTDB (in MIDAS2): Contains 258,406 genomes clustered into 45,555 bacterial and 2,339 archaeal species.

## Example analysis
The function includes filtering of species based on clinical / environmental variables interested, calling of consensus sequencing, constructing the tree, and functional annotations. The below example shows the analysis of the subset of `PRJEB9584`, sequenced by HiSeq 2000. Packages including `ggtree`, `ComplexHeatmap`, `simplyfyEnrichment`, `phangorn` are used.

<img src="https://github.com/noriakis/software/blob/main/images/stana_example.png?raw=true" width="800px">


## TODO
