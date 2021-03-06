% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getMDFR.Surveys.Abundance.R
\name{getMDFR.Surveys.Abundance}
\alias{getMDFR.Surveys.Abundance}
\title{Get survey abundance time series from model results from rsimTCSAM model runs as a dataframe}
\usage{
getMDFR.Surveys.Abundance(rsims, category = "index", cast = "y+x",
  verbose = FALSE)
}
\arguments{
\item{rsims}{- single rsimTCSAM.rep object, rsimTCSAM.resLst object, or named list of the latter}

\item{category}{- 'index' is only choice}

\item{cast}{- casting formula for excluding y,x,m,s,z factor levels from an average-at-size across unspecified factors}

\item{verbose}{- flag (T/F) to print debug info}
}
\value{
dataframe in canonical format
}
\description{
Function to get survey abundance time series from model results from rsimTCSAM model runs as a dataframe.
}
\details{
Extracts the estimated survey abundance time series.
}

