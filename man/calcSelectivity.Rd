% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SelectivityFunctions.R
\name{calcSelectivity}
\alias{calcSelectivity}
\title{Calculate values for a selectivity curve}
\usage{
calcSelectivity(type, z, params)
}
\arguments{
\item{type}{- the type of selectivity function to calculate}

\item{z}{- vector of values at which to calculate the function}

\item{params}{- the selectivity function parameters, as a vector}
}
\value{
vector matching size of z, with names given by elements of z
}

