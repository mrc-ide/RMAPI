% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_map}
\alias{plot_map}
\title{Plot hex map of RMAPI output}
\usage{
plot_map(
  proj,
  variable = NULL,
  col_scale = viridisLite::magma(100),
  barrier_list = list()
)
}
\arguments{
\item{proj}{object of class \code{rmapi_project}.}

\item{variable}{which element of the project output to use as map colours.}

\item{col_scale}{the colour scale to use.}

\item{barrier_list}{optional list of polygon coordinates that are added to
plot.}
}
\description{
Plot hex map of RMAPI output.
}
