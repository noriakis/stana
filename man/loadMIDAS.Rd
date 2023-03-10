% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkProfile.R
\name{loadMIDAS}
\alias{loadMIDAS}
\title{loadMIDAS}
\usage{
loadMIDAS(
  midas_merge_dir,
  cl,
  filtType = "group",
  candSp = NULL,
  filtNum = 2,
  filtPer = 0.8,
  geneType = "copynum"
)
}
\arguments{
\item{midas_merge_dir}{output directory of merge_midas.py}

\item{cl}{named list of sample IDs}

\item{filtType}{"whole" or "group"}

\item{candSp}{candidate species ID}

\item{filtNum}{The species with number above this threshold
for each category is returned}

\item{filtPer}{filter by percentage}

\item{geneType}{"presabs" or "copynum"}
}
\description{
Assess and store profile for species and return filtered species 
based on the number of samples for each category or whole population.
For MIDAS only.
}
