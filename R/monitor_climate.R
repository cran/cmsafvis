#' A 'cmsaf' extension for creating various climate plots.
#'
#' This plotting routine generates graphical output of the evolution of
#' the given variable within the given time range and area. The intended
#' application is for daily accumulated data, such as sunshine duration.
#' Dependent on the output format a PNG or MP4 is created.
#'
#' You can pass a YAML config file and/or specify the arguments directly.
#' Argument prioritization is done in the following way:
#' Direct argument > config argument > default argument.
#' Thus, if you pass a existing config file but also want to modify a specific argument you can do that easily.
#'
#' @param plot_type Specifies the type of the plot
#' ('absolute_map', 'anomaly_map', 'climatology_map', 'fieldmean_plot', or 'fieldmean_and_anomaly_map').
#' @param config Path to YAML config file (character). The config file does not have
#' to specify all arguments. Each addressed argument has to be formatted according to the example config file: (#TODO: LINK EXAMPLE CONFIG FILE!).
#' @param variable Name of variable in infile (\code{NULL} or character). If \code{NULL} then the first variable from the infile is taken.
#' @param accumulate Whether the input file should be accumulated (logical).
#' @param infile Path to NetCDF file (\code{NULL} or character). If \code{NULL} then it needs to be specified in the config file.
#' @param temp_dir Path to temporary working directory (character).
#' @param out_dir Path to output directory (character).
#' @param climate_dir Path to directory in which climatology is computed or contained (\code{NULL} or character). If \code{NULL} then the temp_dir directory is taken.
#' @param climate_year_start Start year of climatology (integer).
#' @param climate_year_end End year of climatology (integer).
#' @param show_extreme_climate_years Whether the minimum and maximum of the climate years should be titled in the fieldmean plot (\code{NULL} or logical).
#' This is usually only of interest when plotting accumulated data. If the default \code{NULL} is chosen, then it will be set to the value of \code{accumulate}.
#' @param climatology_until_eoy Plot the climatology and fieldmeans until the end of year (logical). Only affects fieldmean plots analyzed from January 1st.
#' @param start_date Start date in format of 'YYYY-MM-DD' (\code{NULL} or character). If \code{NULL} then the first date of the infile is used.
#' @param end_date End date in format of 'YYYY-MM-DD' (\code{NULL} or character). If \code{NULL} then the last date of the infile is used.
#' @param country_code Either a country code in iso3c format or from the following: 'AFR' for Africa, 'EUR' for Europe, 'TOT' for the total disc,
#' or 'S_A' for an arbitrary region selection (character). If a country is passed the data from within this country is extracted,
#' else a rectangular box is visualized. Directly provided latitude and longitude ranges will be ignored in case of 'AFR', 'EUR' or 'TOT'.
#' @param lon_min Longitude of lower left corner (\code{NULL} or numeric). If \code{NULL} then the smallest longitude of the infile is used.
#' @param lon_max Longitude of upper right left corner (\code{NULL} or numeric). If \code{NULL} then the largest longitude of the infile is used.
#' @param lat_min Latitude of lower left corner (\code{NULL} or numeric). If \code{NULL} then the smallest latitude of the infile is used.
#' @param lat_max Latitude of upper right corner (\code{NULL} or numeric). If \code{NULL} then the largest latitude of the infile is used.
#' @param outfile_name Filename of the PNG or MP4 outfile (\code{NULL} or character).
#' If \code{NULL} then a name is computed from the current configuration.
#' Please match the file ending according to the output_format.
#' @param output_format Specification of output format (either 'graphic' for PNG or 'animation' for MP4).
#' @param animation_pace Pace of the animation in seconds (positive numeric). This only has an effect if \code{output_format == 'animation'}.
#' @param freeze_animation If TRUE then the animation will freeze at the last frame (logical).
#' @param min_value Lower values than this are ignored (\code{NULL} or numeric). If \code{NULL}, no values are ignored.
#' @param max_value Larger values than this are ignored (\code{NULL} or numeric). If \code{NULL}, no values are ignored.
#' @param nbreaks Number of color breaks (\code{NULL} or positive integer). A value will be computed if \code{NULL} is passed.
#' @param language Language used for title, legend, etc. in plots (either 'eng' for English or 'deu' for German).
#' @param keep_files A flag indicating whether all files created in the process of obtaining the output file should be kept (logical).
#' If false, all intermediate results are deleted, otherwise all are kept. Keeping these files could improve performance in further function calls.
#' @param states Whether to crop/plot administration level of states (logical).
#' @param attach Whether to temporarily merge the infile to an already existing one. (logical).
#' @param infile_attach File to attach the infile to. When 'auto', a suitable file will be searched in out_dir.
#' If attach is false, this will be ignored(character).
#' @param verbose Whether to display progress messages (logical).
#'
#' @export
#' @importFrom assertthat assert_that is.dir is.readable is.string
#' @examples
#'## Create an example NetCDF file with a similar structure as used by CM
#'## SAF. The file is created with the ncdf4 package.  Alternatively
#'## example data can be freely downloaded here: <https://wui.cmsaf.eu/>
#'
#'library(ncdf4)
#'
#'## create some (non-realistic) example data
#'
#'lon <- seq(5, 15, 0.5)
#'lat <- seq(45, 55, 0.5)
#'time <- seq(as.Date("2010-01-01"), as.Date("2011-12-31"), "days")
#'origin <- as.Date("1983-01-01 00:00:00")
#'time <- as.numeric(difftime(time, origin, units = "hour"))
#'data <- array(250:350, dim = c(21, 21, 2 * 365))
#'
#'## create example NetCDF
#'infile <- tempfile("input", fileext = ".nc")
#'
#'x <- ncdim_def(name = "lon", units = "degrees_east", vals = lon)
#'y <- ncdim_def(name = "lat", units = "degrees_north", vals = lat)
#'t <- ncdim_def(name = "time", units = "hours since 1983-01-01 00:00:00",
#'  vals = time, unlim = TRUE)
#'var1 <- ncvar_def("SDU", "W m-2", list(x, y, t), -1, prec = "short")
#'vars <- list(var1)
#'ncnew <- nc_create(infile, vars)
#'ncvar_put(ncnew, var1, data)
#'ncatt_put(ncnew, "lon", "standard_name", "longitude", prec = "text")
#'ncatt_put(ncnew, "lat", "standard_name", "latitude", prec = "text")
#'nc_close(ncnew)
#'
#'## this will save 'output.png' in temp directory
#'cmsafvis::monitor_climate(
#'  plot_type = "anomaly_map",
#'  accumulate = TRUE,
#'  infile = infile,
#'  out_dir = tempdir(),
#'  climate_dir = tempdir(),
#'  climate_year_start = 2010,
#'  climate_year_end = 2011,
#'  end_date = "2010-01-05",
#'  outfile_name = "output",
#'  output_format = "graphic",
#'  keep_files = FALSE,
#'  verbose = FALSE
#')
monitor_climate <- function(plot_type = "absolute_map",
                            config = NULL,
                            variable = NULL,
                            accumulate = FALSE,
                            infile = NULL,
                            temp_dir = tempdir(),
                            out_dir = getwd(),
                            climate_dir = NULL,
                            climate_year_start = 1983,
                            climate_year_end = 2018,
                            show_extreme_climate_years = NULL,
                            climatology_until_eoy = FALSE,
                            start_date = NULL,
                            end_date = NULL,
                            country_code = "S_A",
                            lon_min = NULL,
                            lon_max = NULL,
                            lat_min = NULL,
                            lat_max = NULL,
                            outfile_name = NULL,
                            output_format = "animation",
                            animation_pace = 0.07,
                            freeze_animation = FALSE,
                            min_value = NULL,
                            max_value = NULL,
                            nbreaks = NULL,
                            language = "eng",
                            keep_files = TRUE,
                            states = FALSE,
                            attach = FALSE,
                            infile_attach = "auto",
                            verbose = TRUE
) {

  #### Deal with config ####
  if (is.null(config)) {
    configRead <- FALSE
  } else {
    assert_that(is.string(config))
    assert_that(file.exists(normalizePath(config, mustWork = FALSE)))
    assert_that(!is.dir(config))
    assert_that(is.readable(config))

    configParams <- yaml::read_yaml(config)

    configRead <- TRUE
  }

  #### plot_type ####
  if (missing(plot_type)) {
    if (configRead) {
      plot_type <- tryCatch(
        configParams$plot_type,
        error = function(cond) {
          stop(paste0(
            "Can't read 'plot_type' from config file: ",
            normalizePath(config)
          ))
        }
      )
    } else {
      plot_type <- "absolute_map"
    }
  }
  assert_that(is.string(plot_type))
  assert_that(
    plot_type %in% c(
      "absolute_map",
      "anomaly_map",
      "climatology_map",
      "fieldmean_plot",
      "fieldmean_and_anomaly_map"
    ), msg = paste0("Unknown 'plot_type': ", plot_type, ". See ?monitor_climate for details.")
  )

  #### Actual function call depending on plot_type ####
  arguments_necessary <- methods::formalArgs(plot_type)
  arguments <- as.list(match.call())
  arguments[[1]] <- NULL
  arguments <- arguments[stats::na.omit(match(arguments_necessary, names(arguments)))]
  # We need `getExportedValue` here for the case where
  # `cmsafvis::monitor_climate()` is used without previously calling
  # `library("cmsafvis")`.
  do.call(getExportedValue("cmsafvis", plot_type), arguments, envir = parent.frame())
}