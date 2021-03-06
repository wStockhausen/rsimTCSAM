% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseMC.OutputInfo.R
\name{parseMC.OutputInfo}
\alias{parseMC.OutputInfo}
\title{Parse output info for a fishery catch type or survey from a model configuration file.}
\usage{
parseMC.OutputInfo(rsp, i)
}
\arguments{
\item{rsp}{- parsed text list from model configuration file}

\item{i}{- index to start of output info section in rsp}
}
\value{
list with elements
i - index starting next section
output - list object with output info for abundance, biomass and sizecomps data types
}
\description{
Function to parse output info for a fishery catch type or survey from a model configuration file.
}

