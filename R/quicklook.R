#'Create a quicklook of NetCDF data
#'
#'The function creates a plot of the variables in NetCDF file(s) specified in the config file.
#'Only NetCDF files that conform to the [CMSAF naming convention](https://www.cmsaf.eu/EN/Products/NamingConvention/Naming_Convention_node.html) are supported.
#'
#'@param filelist list of NetCDF file to create plots from (character).
#'@param config filename of configuration file. This may include the directory
#'  (character).
#'@param outpath directory in which to save the output files.
#'  (character).
#'@param jpeg_quality jpeg quality for the image in percent, see [grDevices::jpeg()][grDevices::png()]
#'@param dpi resolution of the image in dots per inch, see [grDevices::jpeg()][grDevices::png()]
#'@param iwidth width of the resulting image in pixels, see [grDevices::jpeg()][grDevices::png()]
#'@param logo logical; should the cmsaf logo be added to the plot?
#'@param copyright logical; should the copyright text be added to the plot?
#'@param bluemarble logical; should the data be plotted onto a NASA bluemarble (only available for MSG/Seviri based data)?
#'   Due to data size this option is not available for the cmsafvis package on CRAN. Please have a look at
#'   our website https://www.cmsaf.eu/R_toolbox
#'
#'@return A jpeg file with the same name as the original NetCDF file.
#'@export
#'@importFrom assertthat assert_that is.count is.flag is.readable is.writeable

quicklook <- function(config,
                      filelist,
                      outpath = getwd(),
                      jpeg_quality = 75,
                      dpi = 72,
                      iwidth = 720,
                      logo = TRUE,
                      copyright = TRUE,
                      bluemarble = FALSE) {
  # Make sure that any user settings are reset when the function exits
  # This is a requirement by CRAN
  oldpar <- graphics::par(no.readonly = TRUE)
  # Warning: In graphics::par(oldpar) : par(new) ohne Plot aufgerufen
  on.exit(suppressWarnings(graphics::par(oldpar)))

  ### check parameters ###
  assert_that(is.string(config))
  assert_that(file.exists(normalizePath(config, mustWork = FALSE)))
  assert_that(!is.dir(config))
  assert_that(is.readable(config))

  assert_that(all(file.exists(filelist)))
  for (file_ in filelist) {
    assert_that(is.readable(file_))
  }

  assert_that(is.dir(outpath))
  assert_that(is.writeable(outpath))

  assert_that(is.count(jpeg_quality))
  assert_that(0 <= jpeg_quality && jpeg_quality <= 100)

  assert_that(is.count(dpi))
  assert_that(is.count(iwidth))
  assert_that(is.flag(logo))
  assert_that(is.flag(copyright))
  assert_that(is.flag(bluemarble))

  ### Build colorpalettes ###

  palettes <- GetPaletteConfig(gui = TRUE)
  names(palettes) <- tolower(names(palettes))
  names(palettes)[names(palettes) == "typ"] <- "type"

  # add more color schemes
  new_row <- data.frame("more", NA, NA, NA, NA, NA, NA, NA, NA, NA, 1)
  names(new_row) <- names(palettes)
  palettes <- rbind(palettes, new_row)
  rownames(palettes)[75] <- "tim.colors"

  palettes <- rbind(palettes, new_row)
  rownames(palettes)[76] <- "sunny"
  palettes <- rbind(palettes, new_row)
  rownames(palettes)[77] <- "cloud_mask1"
  palettes <- rbind(palettes, new_row)
  rownames(palettes)[78] <- "cloud_mask2"

  cloud_mask1 <- c("black", "transparent", "gray60", "white")
  cloud_mask2 <- c("black", "transparent", "gray60", "white", "pink")

  ### Read and format logo ###

  if (logo) {
    lf_black <- ifelse(logo == "color","CMSAF_NoName_Colour_crop.png","CMSAF_NoName_Black.png")

    logo_cmsaf_path_black <- system.file(
      "extdata",
      lf_black,
      package = "cmsafvis",
      mustWork = TRUE
    )

    # Size and location of the logo
    logo.scale_black <- 0.3

    logo.x <- 0
    logo.y <- 0

    logo_cmsaf_black <- png::readPNG(logo_cmsaf_path_black)

    dims_black <- dim(logo_cmsaf_black)[1:2]

    logo.height_black <- logo.scale_black * iwidth * dims_black[1] / dims_black[2]

    if (logo.scale_black * iwidth * dims_black[1] / dims_black[2] < 45) {
      logo.height_black <- 45
      logo.scale_black <- logo.height_black / dims_black[1] * dims_black[2] / iwidth
    }

    # Prepare color logo
    lf_color <- "CMSAF_NoName_Colour_crop.png"

    logo_cmsaf_path_color <- system.file(
      "extdata",
      lf_color,
      package = "cmsafvis",
      mustWork = TRUE
    )

    # Size and location of the logo
    logo.scale_color <- 0.3

    logo.x <- 0
    logo.y <- 0

    logo_cmsaf_color <- png::readPNG(logo_cmsaf_path_color)

    dims_color <- dim(logo_cmsaf_color)[1:2]

    logo.height_color <- logo.scale_color * iwidth * dims_color[1] / dims_color[2]

    if (logo.scale_color * iwidth * dims_color[1] / dims_color[2] < 45) {
      logo.height_color <- 45
      logo.scale_color <- logo.height_color / dims_color[1] * dims_color[2] / iwidth
    }

    text.x <- 0.99
    text.y <- 0.01
  }

  ### Read config file ###

  configParams <- yaml::read_yaml(config)
  ref_file <- filelist[1]
  file_info <- get_file_info(ref_file)

  varnames <- c()
  units <- c()
  stacks <- c()
  plot_lim <- c()
  col_from_config <- c()
  legends <- c()
  logos <- c()

  if (bluemarble && !file_info$grid == "Satellite projection MSG/Seviri") {
    stop("Bluemarble plotting is only available for CLAAS data on MSG grid.")
  }

  vars <- names(configParams[[file_info$product_type]][[file_info$id]])[names(configParams[[file_info$product_type]][[file_info$id]]) != "Dataset"]
  nvars <- length(vars)

  # no plot variables found
  assert_that(nvars > 0)  # TODO Improve error message so that user knows how to fix the problem.
  is_multiplot <- nvars > 1

  dataset_name <- configParams[[file_info$product_type]][[file_info$id]]$Dataset

  for (i in seq_along(vars)) {
    limits <- configParams[[file_info$product_type]][[file_info$id]][[vars[i]]]$limits
    plot_lim <- rbind(plot_lim, c(limits$min, limits$max))
    legends <- c(legends, configParams[[file_info$product_type]][[file_info$id]][[vars[i]]]$legend)
    col_from_config <- c(col_from_config, configParams[[file_info$product_type]][[file_info$id]][[vars[i]]]$colorscale)
    if (logo) logos <- c(logos, configParams[[file_info$product_type]][[file_info$id]][[vars[i]]]$logo)
  }


  ### Read infiles ###

  nc <- ncdf4::nc_open(ref_file)
  vars <- names(nc$var)[toupper(names(nc$var)) %in% vars]
  nvars <- length(vars)

  # no plot variables found
  assert_that(nvars > 0)  # TODO Improve error message so that user knows how to fix the problem.

  for (k in seq_along(vars)) {
    varnames <- c(varnames, ncdf4::ncatt_get(nc, vars[k], "long_name")$value)
    units <- c(units, ncdf4::ncatt_get(nc, vars[k], "units")$value)
    stacks <- c(stacks, raster::stack(filelist, quick = TRUE, varname = vars[k]))
  }
  if (ncdf4::ncatt_get(nc, 0, "CMSAF_proj4_params")$hasatt) {
    nc_crs <- raster::crs(ncdf4::ncatt_get(nc, 0, "CMSAF_proj4_params")$value)
    for (l in seq_along(stacks)) {
      raster::crs(stacks[[l]]) <- nc_crs
    }
  }
  if ("lon" %in% names(nc$dim)) {
    lon_min <- min(ncdf4::ncvar_get(nc, "lon"), na.rm = TRUE)
    lon_max <- max(ncdf4::ncvar_get(nc, "lon"), na.rm = TRUE)
    lat_min <- min(ncdf4::ncvar_get(nc, "lat"), na.rm = TRUE)
    lat_max <- max(ncdf4::ncvar_get(nc, "lat"), na.rm = TRUE)
  } else if (ncdf4::ncatt_get(nc, 0, "geospatial_lon_max")$hasatt) {
    lon_min <- ncdf4::ncatt_get(nc, 0, "geospatial_lon_min")$value
    lon_max <- ncdf4::ncatt_get(nc, 0, "geospatial_lon_max")$value
    lat_min <- ncdf4::ncatt_get(nc, 0, "geospatial_lat_min")$value
    lat_max <- ncdf4::ncatt_get(nc, 0, "geospatial_lat_max")$value
  } else {
    stop("unable to get a lon/lat reference")
  }

  for (l in seq_along(stacks)) {
    raster::extent(stacks[[l]]) <- c(lon_min, lon_max, lat_min, lat_max)
  }
  ncdf4::nc_close(nc)

  ### Set aspect ###

  lon_range <- lon_max - lon_min
  lat_range <- lat_max - lat_min

  aspect <- lon_range/lat_range

  iheight <- round(iwidth/aspect)

  if (logo) {
    AR_color <- dims_color[1] / dims_color[2] * iwidth / iheight
    AR_black <- dims_black[1] / dims_black[2] * iwidth / iheight
  }

  ### Plot ###

  for (i in 1:raster::nlayers(stacks[[1]])) {
    # filename and timestamp for title
    filename <- unlist(strsplit(basename(filelist[i]), "\\."))
    outfile <- file.path(outpath, paste0(filename[1], ".jpg"))
    fi <- get_file_info(filename[1])
    file_time <- fi$date_time
    if (file_info$time_interval == "instantaneous")
      file_time <- format(file_time, "%Y-%m-%d %R")

    grDevices::jpeg(outfile,
                    quality = jpeg_quality,
                    width = iwidth,
                    height = iheight,
                    res = dpi,
                    pointsize = round(12 * (1/(dpi/72)))
                    )

    # Parameters for multiple plots

    if (is_multiplot) {
      ncols <- ceiling(sqrt(nvars))
      nrows <- ceiling(nvars/ncols)
      graphics::par(mfrow = c(nrows, ncols))
    }

    graphics::par(mar = c(5, 4, 6, 4) + 0.1)

    for (j in seq_along(vars)) {

      # Set color palette
      if (col_from_config[[1]] == "clouds") {
        stacks[[j]][[i]][is.na(stacks[[j]][[i]])] <- 0
        if (raster::maxValue(stacks[[j]][[i]]) == 3) {
          col <- cloud_mask1
        } else {
          col <- cloud_mask2
        }
        plot_lim[j,] <- range(raster::values(stacks[[j]][[i]]))
      } else {
        col <- getColors(col_from_config[[j]], palettes, 32, FALSE)
      }
      # bluemarble plot
      if (bluemarble) {
		stop("Bluemarble plotting is not available. See https://www.cmsaf.eu/R_toolbox")
        # raster::extent(stacks[[j]]) <- c(-1, 1, -1, 1)
        # fields::quilt.plot(
          # # This is generated in data-raw/generate_internal_data.R
          # blue_marble$projection$x,
          # blue_marble$projection$y,
          # blue_marble$data_values,
          # xlim = c(-1, 1),
          # ylim = c(-1, 1),
          # nx = blue_marble$n_lon_unique / blue_marble$xf,
          # ny = blue_marble$n_lat_unique / blue_marble$yf,
          # xlab = " ",
          # ylab = " ",
          # main = "",
          # col = blue_marble$colors,
          # add.legend = FALSE,
          # axes = FALSE
        # )
        # raster::image(stacks[[j]], y = i,
                      # main = "",
                      # axes = FALSE,
                      # xlab = "",
                      # ylab = "",
                      # col = col,
                      # add = TRUE)

      } else {
        # borderline plots for scale
        if (file_info$grid == "Satellite projection MSG/Seviri") {
          stop("Bluemarble plotting is not available. See https://www.cmsaf.eu/R_toolbox")
		  # graphics::par(pty = "s")
          # raster::extent(stacks[[j]]) <- c(-1, 1, -1, 1)
          # fields::quilt.plot(
            # # This is generated in data-raw/generate_internal_data.R
            # blue_marble$projection$x,
            # blue_marble$projection$y,
            # blue_marble$data_values,
            # xlim = c(-1, 1),
            # ylim = c(-1, 1),
            # nx = blue_marble$n_lon_unique / blue_marble$xf,
            # ny = blue_marble$n_lat_unique / blue_marble$yf,
            # xlab = " ",
            # ylab = " ",
            # col = "gray",
            # add.legend = FALSE,
            # axes = FALSE
          # )
        } else {
          graphics::image(lon_min:(lon_max*1.2),
                          lat_min:lat_max,
                          outer(lon_min:(lon_max*1.2),lat_min:lat_max,"+"),
                          main = "",
                          xlim = c(lon_min, lon_max),
                          ylim = c(lat_min, lat_max),
                          xlab = " ",
                          ylab = " ",
                          col = "gray",
                          axes = FALSE
          )
          }

        # plot image
        raster::image(stacks[[j]], y = i,
                      main = "",
                      xlim = c(lon_min, lon_max),
                      ylim = c(lat_min, lat_max),
                      axes = FALSE,
                      xlab = "",
                      ylab = "",
                      zlim = plot_lim[j,],
                      col = col,
                      colNA = "gray",
                      add = TRUE
        )

        # borderline plot
        if (file_info$grid == "Satellite projection MSG/Seviri") {
          raster::extent(stacks[[j]]) <- c(-1, 1, -1, 1)
          suppressWarnings(
            maps::map("world", projection = "orthographic", interior = FALSE, orientation = c(0,0,0), add = TRUE)
          )
        } else {
          maps::map("world", interior = FALSE, xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), add = TRUE)
        }
}
        # plot logo and copyright text
        if (logo || copyright) {
          graphics::par(usr = c(0, 1, 0, 1))

          if (logos[j] == "color") {
            logo_cmsaf <- logo_cmsaf_color
            AR <- AR_color
            logo.scale <- logo.scale_color
            logo.height <- logo.height_color
          } else {
            logo_cmsaf <- logo_cmsaf_black
            AR <- AR_black
            logo.scale <- logo.scale_black
            logo.height <- logo.height_black
          }
        }
        if (logo) {
          graphics::rasterImage(array(0.75, dim = dim(logo_cmsaf)),
                                logo.x + 0.01,
                                logo.y + 0.01,
                                logo.x + logo.scale,
                                logo.y + (AR * logo.scale),
                                interpolate = TRUE,
                                bg = "white")

          graphics::rasterImage(logo_cmsaf,
                                logo.x + 0.01,
                                logo.y + 0.01,
                                logo.x + logo.scale,
                                logo.y + (AR * logo.scale),
                                interpolate = TRUE
          )
        }
        if (copyright) {
          txt <- paste0("\u00a9 EUMETSAT, ", format.Date(Sys.Date(), "%Y"))
          cr.scale <- (logo.scale - 0.02)/graphics::strwidth(txt)
          dims_text <- c(round(dim(logo_cmsaf)[1]/2), dim(logo_cmsaf)[2], dim(logo_cmsaf)[3])
          graphics::rasterImage(array(0.75, dim = dim(logo_cmsaf)),
                                text.x - logo.scale,
                                text.y,
                                text.x,
                                text.y + ((AR * logo.scale)/2),
                                interpolate = TRUE)

          graphics::text(text.x, text.y + 0.01, txt,
                         cex = cr.scale,
                         adj = c(1,0))
        }

        # plot legend
        if (legends[j]) {
          raster::plot(stacks[[j]], y = i,
                       main = "",
                       axes = FALSE,
                       xlab = "",
                       ylab = "",
                       zlim = plot_lim[j,],
                       legend.only = TRUE,
                       legend.shrink = 1,
                       legend.width = 1,
                       legend.mar = 4,
                       legend.lab = units[j],
                       col = col,
                       add = TRUE)
        }


      # figure title
      if (is_multiplot) graphics::mtext(varnames[j], line = 1, cex = 1.6)

    }

    # main title
    if (is_multiplot) {
      graphics::mtext(paste(dataset_name, file_time, sep = ", "),
                      side = 3,
                      line = -2,
                      outer = TRUE,
                      cex = 2,
                      font = 2
      )
    } else {
      graphics::mtext(paste(varnames[1], dataset_name, file_time, sep = ", "),
                      side = 3,
                      line = -2,
                      outer = TRUE,
                      cex = 2,
                      font = 2
      )
    }
    grDevices::dev.off()
  }
}
