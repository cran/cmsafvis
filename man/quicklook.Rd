% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quicklook.R
\name{quicklook}
\alias{quicklook}
\title{Create a quicklook of NetCDF data}
\usage{
quicklook(
  config,
  filelist,
  outpath = getwd(),
  jpeg_quality = 75,
  dpi = 72,
  iwidth = 720,
  logo = TRUE,
  copyright = TRUE,
  bluemarble = FALSE
)
}
\arguments{
\item{config}{filename of configuration file. This may include the directory
(character).}

\item{filelist}{list of NetCDF file to create plots from (character).}

\item{outpath}{directory in which to save the output files.
(character).}

\item{jpeg_quality}{jpeg quality for the image in percent, see \link[grDevices:png]{grDevices::jpeg()}}

\item{dpi}{resolution of the image in dots per inch, see \link[grDevices:png]{grDevices::jpeg()}}

\item{iwidth}{width of the resulting image in pixels, see \link[grDevices:png]{grDevices::jpeg()}}

\item{logo}{logical; should the cmsaf logo be added to the plot?}

\item{copyright}{logical; should the copyright text be added to the plot?}

\item{bluemarble}{logical; should the data be plotted onto a NASA bluemarble (only available for MSG/Seviri based data)?
Due to data size this option is not available for the cmsafvis package on CRAN. Please have a look at
our website https://www.cmsaf.eu/R_toolbox}
}
\value{
A jpeg file with the same name as the original NetCDF file.
}
\description{
The function creates a plot of the variables in NetCDF file(s) specified in the config file.
Only NetCDF files that conform to the \href{https://www.cmsaf.eu/EN/Products/NamingConvention/Naming_Convention_node.html}{CMSAF naming convention} are supported.
}
