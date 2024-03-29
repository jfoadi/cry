% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_SHELX.R
\name{plot_SHELX}
\alias{plot_SHELX}
\title{Plot SHELXC log files}
\usage{
plot_SHELX(filename, filename_e, var, type, title_chart)
}
\arguments{
\item{filename}{A data frame in output from \code{{read_SHELX_log}}.}

\item{filename_e}{A data frame with the inverted hand from shelxe}

\item{var}{the variable to be plotted vs the resolution}

\item{type}{indicate the type of file, possible value are "shelxc",
"shelxd" and "shelxe".}

\item{title_chart}{title of the chart.}
}
\value{
A graphical object from ggplot2 class that contains
the solution founded by SHELX log file.
}
\description{
Plot SHELXC log files
}
\examples{
datadir <- system.file("extdata",package="cry")
## SHELXC
shelxc_log <- file.path(datadir,"shelxc.log")
shelxc <- read_SHELX_log(shelxc_log)
plot_shelxc <- plot_SHELX(filename = shelxc, var = shelxc$I_sig,
type = "shelxc", title_chart = "SHELXC")
plot_shelxc
## SHELXD
shelxd_log <- file.path(datadir,"shelxd.log")
shelxd <- read_SHELX_log(shelxd_log)
plot_shelxd <- plot_SHELX(filename = shelxd, type = "shelxd",
title_chart = "SHELXD")
plot_shelxd
## SHELXE
filename_i <- file.path(datadir,"shelxe_i.log")
shelxe_i <- read_SHELX_log(filename_i)
filename_o <- file.path(datadir,"shelxe_o.log")
shelxe_o <- read_SHELX_log(filename_o)
plot_shelxe <- plot_SHELX(filename = shelxe_i,
filename_e = shelxe_o, type = "shelxe", title_chart = "SHELXE")
plot_shelxe

}
