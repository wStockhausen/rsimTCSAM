% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel.R
\name{runModel}
\alias{runModel}
\title{Run the simulation model}
\usage{
runModel(mc, mp, showPlot = TRUE)
}
\arguments{
\item{mc}{- model configuration object}

\item{mp}{- model processes object}

\item{showPlot}{- flag to show plots}
}
\value{
List consisting of:
iN_xmsz - initial population abundance by sex/maturity/shell condition/size
P_list - list of population time series (see calcNatZ)
F_list - list of fisheries results (see calcNatZ.Fisheries)
S_list - list of survey results (see calcNatZ.Surveys)
}
\description{
Function to run the simulation model
}

