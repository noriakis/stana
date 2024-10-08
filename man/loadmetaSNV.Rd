% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadmetaSNV.R
\name{loadmetaSNV}
\alias{loadmetaSNV}
\title{loadmetaSNV}
\usage{
loadmetaSNV(metasnv_out_dir, cl = NULL, just_species = FALSE, candSp = NULL)
}
\arguments{
\item{metasnv_out_dir}{output directory of merge_midas.py}

\item{cl}{named list of sample IDs}

\item{just_species}{just return species id}

\item{candSp}{candidate species ID}
}
\value{
stana object
}
\description{
Assess and store profile for species for metaSNV.
}
