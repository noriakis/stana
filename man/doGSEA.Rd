% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doGSEA.R
\name{doGSEA}
\alias{doGSEA}
\title{doGSEA}
\usage{
doGSEA(
  stana,
  candSp = NULL,
  cl = NULL,
  eps = 0.01,
  how = sum,
  zeroPerc = 0,
  rankMethod = "modt",
  target = "pathway",
  bg_filter = TRUE,
  background = NULL,
  gseaArgs = list()
)
}
\arguments{
\item{stana}{stana object}

\item{candSp}{candidate species}

\item{cl}{if NULL, grouping in stana is used}

\item{eps}{pseudovalue added if log is taken}

\item{how}{how to combine multiple variables, default to sum}

\item{zeroPerc}{genes >= the percentage of count zero sample will be excluded.
Default to zero, not recommended in GSEA}

\item{rankMethod}{how to rank genes}

\item{target}{pathway or module}

\item{bg_filter}{filter the background for those in table}

\item{background}{background TERM2GENE data}

\item{gseaArgs}{list of argument passed to GSEA or fgsea function}
}
\value{
stana object with GSEA results from clusterProfiler
}
\description{
Based on KEGG database, GSEA will be performed by the function.
clusterProfiler::GSEA is used.
By default this uses moderated t-statistics for ranking of the genes.
}
