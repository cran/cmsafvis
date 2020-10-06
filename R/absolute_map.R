#' A 'cmsaf' extension for creating absolute valued plots.
#'
#' This plotting routine generates graphical output for
#' the given variable within the given time range and area.
#' Dependent on the output format a PNG or MP4 is created.
#'
#' You can pass a YAML config file and/or specify the arguments directly.
#' Argument prioritization is done in the following way:
#' Direct argument > config argument > default argument.
#' Thus, if you pass a existing config file but also want to modify a specific
#' argument you can do that easily.
#'
#' @inheritParams monitor_climate
#'
#' @export
#' @importFrom assertthat assert_that is.date is.dir is.flag is.number is.readable is.string is.writeable
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
#'time <- seq(as.Date("2010-01-01"), as.Date("2010-01-01"), "days")
#'origin <- as.Date("1983-01-01 00:00:00")
#'time <- as.numeric(difftime(time, origin, units = "hour"))
#'data <- array(250:350, dim = c(21, 21, 1))
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
#'cmsafvis::absolute_map(
#'  infile = infile,
#'  out_dir = tempdir(),
#'  outfile_name = "output",
#'  output_format = "graphic",
#'  keep_files = FALSE,
#'  verbose = FALSE
#')
absolute_map <- function(config = NULL,
                         variable = NULL,
                         accumulate = FALSE,
                         infile = NULL,
                         temp_dir = tempdir(),
                         out_dir = getwd(),
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
                         verbose = TRUE) {
  # Call central argument parser
  arguments_necessary <- methods::formalArgs(parse_arguments)
  arguments <- as.list(match.call())
  arguments[[1]] <- "absolute_map"
  names(arguments)[1] <- "plot_type"
  arguments <- arguments[stats::na.omit(match(arguments_necessary, names(arguments)))]
  parsedArguments <- do.call(parse_arguments, arguments, envir = parent.frame())

  # Use parsed arguments
  variable <- parsedArguments$variable
  accumulate <- parsedArguments$accumulate
  temp_dir <- parsedArguments$temp_dir
  country_code <- parsedArguments$country_code
  out_dir <- parsedArguments$out_dir
  lon_min <- parsedArguments$lon_min
  lon_max <- parsedArguments$lon_max
  lat_min <- parsedArguments$lat_min
  lat_max <- parsedArguments$lat_max
  infile <- parsedArguments$infile
  start_date <- parsedArguments$start_date
  end_date <- parsedArguments$end_date
  language <- parsedArguments$language
  animation_pace <- parsedArguments$animation_pace
  output_format <- parsedArguments$output_format
  min_value <- parsedArguments$min_value
  max_value <- parsedArguments$max_value
  nbreaks <- parsedArguments$nbreaks
  freeze_animation <- parsedArguments$freeze_animation
  outfile_name <- parsedArguments$outfile_name
  keep_files <- parsedArguments$keep_files
  states <- parsedArguments$states
  attach <- parsedArguments$attach
  infile_attach <- parsedArguments$infile_attach
  new_infile <- parsedArguments$new_infile
  verbose <- parsedArguments$verbose

  if (attach) {
    attach_file(variable = variable,
                infile = infile,
                infile_attach = infile_attach,
                new_infile = new_infile,
                temp_dir = temp_dir,
                verbose = verbose
    )
    infile <- new_infile
  }

  finalInfile <- extractFinalOutfile(
    variable = variable,
    infile = infile,
    start_date = start_date,
    end_date = end_date,
    accumulate = accumulate,
    temp_dir = temp_dir,
    verbose = verbose
  )

  if (is_country(country_code)) {
    mask_file <- create_country_mask(
      infile = finalInfile,
      temp_dir = temp_dir,
      country_code = country_code,
      states = states,
      verbose = verbose
    )

    mask_file_final <- create_country_mask_final(
      mask_infile = mask_file,
      temp_dir = temp_dir,
      country_code = country_code,
      lon_min = lon_min,
      lon_max = lon_max,
      lat_min = lat_min,
      lat_max = lat_max,
      verbose = verbose
    )

    # Remove reusable file if desired
    if (!keep_files & file.exists(mask_file)) {
      file.remove(mask_file)
    }
  }

  var_masked_file <- absolute_map_builder(
    variable = variable,
    temp_dir = temp_dir,
    infile = finalInfile,
    mask_file_final = mask_file_final,
    country_code = country_code,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    start_date = start_date,
    end_date = end_date
  )

  # Remove reusable file if desired
  if (!keep_files) {
    if (is_country(country_code) && file.exists(mask_file_final)) { file.remove(mask_file_final) }
  }

  plot_abs_map(
    variable = variable,
    country_code = country_code,
    infile = var_masked_file,
    out_dir = out_dir,
    start_date = start_date,
    end_date = end_date,
    language = language,
    plot_type = "absolute_map",
    animation_pace = animation_pace,
    output_format = output_format,
    min_value = min_value,
    max_value = max_value,
    nbreaks = nbreaks,
    freeze_animation = freeze_animation,
    outfile_name = outfile_name,
    accumulate = accumulate,
    states = states,
    verbose = verbose
  )
}