% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NMF.R
\name{plotAbundanceWithinSpecies}
\alias{plotAbundanceWithinSpecies}
\title{plotAbundanceWithinSpecies}
\usage{
plotAbundanceWithinSpecies(
  stana,
  species,
  tss = TRUE,
  return_data = FALSE,
  by = "NMF"
)
}
\arguments{
\item{stana}{stana object}

\item{species}{species ID}

\item{tss}{perform total sum scaling}

\item{return_data}{return only the data, not plot}

\item{by}{NMF or coef matrix set to `coefMat` slot}
}
\value{
ggplot object
}
\description{
plot abundances using factor to sample matrix produced by NMF.
}
