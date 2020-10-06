#' A 'cmsaf' extension for creating both, a spatial mean and an anomaly map.
#'
#' This plotting routine generates a graph showing the evolution of
#' the spatial mean of a given variable and the corresponding anomaly
#' map within the given time range and area. The intended application
#' is for daily accumulated data, such as sunshine duration.
#' Dependent on the output format a PNG or MP4 is created.
#'
#' You can pass a YAML config file and/or specify the arguments directly.
#' Argument prioritization is done in the following way:
#' Direct argument > config argument > default argument.
#' Thus, if you pass a existing config file but also want to modify a specific argument you can do that easily.
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
#'cmsafvis::fieldmean_and_anomaly_map(
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
fieldmean_and_anomaly_map <- function(config = NULL,
                                      variable = NULL,
                                      accumulate = FALSE,
                                      infile = NULL,
                                      temp_dir = tempdir(),
                                      out_dir = getwd(),
                                      climate_dir = NULL,
                                      climate_year_start,
                                      climate_year_end,
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
                                      verbose = TRUE) {
  # Call central argument parser
  arguments_necessary <- methods::formalArgs(parse_arguments)
  arguments <- as.list(match.call())
  arguments[[1]] <- "fieldmean_and_anomaly_map"
  names(arguments)[1] <- "plot_type"
  arguments <- arguments[stats::na.omit(match(arguments_necessary, names(arguments)))]
  parsedArguments <- do.call(parse_arguments, arguments, envir = parent.frame())

  # Use parsed arguments
  variable <- parsedArguments$variable
  accumulate <- parsedArguments$accumulate
  temp_dir <- parsedArguments$temp_dir
  infile <- parsedArguments$infile
  climate_dir <- parsedArguments$climate_dir
  climate_year_start <- parsedArguments$climate_year_start
  climate_year_end <- parsedArguments$climate_year_end
  show_extreme_climate_years <- parsedArguments$show_extreme_climate_years
  climatology_until_eoy <- parsedArguments$climatology_until_eoy
  country_code <- parsedArguments$country_code
  out_dir <- parsedArguments$out_dir
  lon_min <- parsedArguments$lon_min
  lon_max <- parsedArguments$lon_max
  lat_min <- parsedArguments$lat_min
  lat_max <- parsedArguments$lat_max
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

  # If user wants to accumulate the file, we do so here.
  finalInfile <- extractFinalOutfile(
    variable = variable,
    infile = infile,
    start_date = start_date,
    end_date = end_date,
    accumulate = accumulate,
    temp_dir = temp_dir,
    verbose = verbose
  )

  climatology_file <- calculate_climatology(
    variable = variable,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    start_date = start_date,
    end_date = end_date,
    country_code = country_code,
    climate_dir = climate_dir,
    infile = infile,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    accumulate = accumulate,
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

    # Apply final mask file to climatology file
    climate_masked_file <- apply_mask_clima(
      variable = variable,
      mask_file_final = mask_file_final,
      climatology_file = climatology_file,
      temp_dir = temp_dir,
      country_code = country_code,
      climate_year_start = climate_year_start,
      climate_year_end = climate_year_end,
      accumulate = accumulate
    )
  } else {
    climate_masked_file <- climatology_file
  }


  var_climatology_accumulated_fldmean <- fieldmean_climate(
    variable = variable,
    temp_dir = temp_dir,
    climate_infile = climate_masked_file,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    country_code = country_code,
    accumulate = accumulate
  )

  infile_current <- fieldmean_current(
    variable = variable,
    temp_dir = temp_dir,
    infile = finalInfile,
    mask_file_final = mask_file_final,
    country_code = country_code,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    start_date = start_date,
    end_date = end_date
  )

  diffclim_accumulated_file <- apply_mask_current(
    variable = variable,
    temp_dir = temp_dir,
    infile = finalInfile,
    mask_file_final = mask_file_final,
    climatology_file = climatology_file,
    country_code = country_code,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    start_date = start_date,
    end_date = end_date
  )

  fieldmean_ensemble(
    variable = variable,
    infile = infile,
    mask_file_final = mask_file_final,
    temp_dir = temp_dir,
    climate_dir = climate_dir,
    country_code = country_code,
    lon_min = lon_min,
    lon_max = lon_max,
    lat_min = lat_min,
    lat_max = lat_max,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    accumulate = accumulate,
    verbose = verbose,
    keep_files = keep_files
  )

  # Remove reusable files if desired
  if (!keep_files) {
    if (file.exists(climate_masked_file)) { file.remove(climate_masked_file) }
    if (file.exists(climatology_file)) { file.remove(climatology_file) }
    if (is_country(country_code) && file.exists(mask_file_final)) { file.remove(mask_file_final) }
  }

  plot_fieldmean_and_map(
    variable = variable,
    country_code = country_code,
    climate_year_start = climate_year_start,
    climate_year_end = climate_year_end,
    show_extreme_climate_years = show_extreme_climate_years,
    climatology_until_eoy = climatology_until_eoy,
    infile = var_climatology_accumulated_fldmean,
    infile2 = infile_current,
    infile_map = diffclim_accumulated_file,
    temp_dir = temp_dir,
    out_dir = out_dir,
    start_date = start_date,
    end_date = end_date,
    language = language,
    animation_pace = animation_pace,
    output_format = output_format,
    min_value = min_value,
    max_value = max_value,
    nbreaks = nbreaks,
    freeze_animation = freeze_animation,
    outfile_name = outfile_name,
    adjustAccumulation = accumulate,
    states = states,
    verbose = verbose
  )
}