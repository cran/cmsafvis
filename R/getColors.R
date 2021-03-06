# Function for getting colors. Either the basic color schemes or from colorspace.
getColors <- function(PAL,
                      palettes,
                      num_brk,
                      reverse) {
  idx <- which(rownames(palettes) == PAL)
  name   <- PAL
  if (PAL == "tim.colors") {
    name <- "fields::tim.colors"
  }
  if (PAL == "bpy") {
    name <- "sp::bpy.colors"
  }
  curPAL <- as.list(palettes[idx, ])
  if (length(idx) == 0) {
    idx <- which(rownames(palettes) == "sunny")
    name   <- "sunny"
    curPAL <- as.list(palettes[idx, ])
  }

  if (curPAL$type == "base") {
    pal <- eval(parse(text = tolower(name)))
  } else if (curPAL$type == "more") {
    # somehow sunny IS used
    sunny <- grDevices::colorRampPalette(c("black",
                                           "#3a0303",
                                           "#640000",
                                           "#981800",
                                           "#ca4b00",
                                           "#fc7f01",
                                           "#ffb234",
                                           "#ffe566",
                                           "#ffff98",
                                           "#ffffcb",
                                           "white"))
    cloud_mask1 <- grDevices::colorRampPalette(c("black", "transparent", "gray60", "white"))
    cloud_mask2 <- grDevices::colorRampPalette(c("black", "transparent", "gray60", "white", "pink"))
    pal <- eval(parse(text = tolower(name)))
  } else {
    curPAL$reverse <- FALSE
    pal <- do.call(GetPalette, curPAL)
  }

  colorbar <- pal(num_brk)

  if (reverse) {
    colorbar <- rev(colorbar)
  }

  return(colorbar)
}
