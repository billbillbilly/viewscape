#' @import sf
#' @importFrom graphics par
#' @importFrom grDevices rgb
#' @importFrom methods new
#' @importFrom stats sd
#' @importFrom sp SpatialPoints
#' @importFrom terra as.lines

#' @noMd
# Convert a value in metres to the native linear unit of a SpatRaster CRS.
# Feet-based CRS identifiers include "ft", "us-ft", "foot", and "feet".
# If the unit cannot be determined, the value is returned unchanged (assumes metres).
m_to_crs_units <- function(value_m, dsm) {
  units <- sf::st_crs(dsm)$units
  if (!is.null(units) && tolower(units) %in% c("ft", "us-ft", "foot", "feet")) {
    return(value_m * 3.28084)
  }
  return(value_m)
}

#' @noMd
radius_viewshed <- function(dsm, r, refraction_factor, viewPt, offset, offset2 = 0, method) {
  resolution <- terra::res(dsm)
  distance <- round(r/resolution[1])
  projection <- terra::crs(dsm, proj = TRUE)
  # create an extent to crop input raster
  subarea <- get_buffer(viewPt[1], viewPt[2], r)
  subdsm <- terra::crop(dsm, terra::ext(subarea))
  dsm <- subdsm
  # setup the view point
  col <- terra::colFromX(dsm, viewPt[1])
  row <- terra::rowFromY(dsm, viewPt[2])
  .ex_vp = terra::extract(dsm,cbind(viewPt[1],viewPt[2]))
  z_viewpoint = .ex_vp[1,ncol(.ex_vp)]+offset
  viewpoint <- matrix(0,1,3)
  viewpoint[1,1] <- col
  viewpoint[1,2] <- row
  viewpoint[1,3] <- z_viewpoint
  # get raster information
  dsm_matrix <- terra::as.matrix(dsm, wide=TRUE)
  # compute viewshed
  if (method == "plane") {
    label_matrix <- reference(viewpoint, dsm_matrix, offset2, distance, refraction_factor)
  } else if (method == "los") {
    label_matrix <- LOS(viewpoint, dsm_matrix, offset2, distance, refraction_factor)
  }

  output <- new("Viewshed",
                viewpoint = c(viewPt, offset),
                viewpos = c(col,row),
                visible = label_matrix,
                resolution = resolution,
                extent = as.vector(sf::st_bbox(dsm)),
                crs = projection)
  return(output)
}

#' @noMd
paral_nix <- function(X, dsm, r, refraction_factor, offset, workers, method){
  results <- pbmcapply::pbmclapply(X = X,
                                   FUN=radius_viewshed,
                                   dsm=dsm,
                                   r=r,
                                   refraction_factor=refraction_factor,
                                   offset=offset,
                                   method = method,
                                   mc.cores=workers)
  return(results)
}

#' @noMd
# H=−∑[(pi)×ln(pi)]
sd_index <- function(p) {
  out <- sum(log(p) * p) * -1
  return(round(out, digits = 3))
}

#' @noMd
# create a buffer based on a given point
get_buffer <- function(x, y, r) {
  pdf <- data.frame(row.names = 1)
  pdf[1,"x"] <- x
  pdf[1,"y"] <- y
  p <- sp::SpatialPoints(pdf)
  p <- sf::st_as_sf(p)
  subarea <- sf::st_buffer(p, r)
  return(subarea)
}

#' @noMd
# get patches
get_patch <- function(viewshed){
  vpt <- filter_invisible(viewshed, FALSE)
  m <- terra::vect(sp::SpatialPoints(vpt))
  terra::crs(m) <- viewshed@crs
  mask_ <- terra::mask(filter_invisible(viewshed, TRUE), m)
  return(mask_)
}

#' @noMd
# get patches parameter
patch_p <- function(m, patchpoly){
  clusters <- terra::patches(m, directions=4)
  ptc <- terra::as.polygons(clusters)
  patchpoly <- terra::mask(patchpoly, ptc)
  ptc_lines <- m %>%
    terra::as.polygons() %>%
    terra::as.lines() %>%
    sf::st_as_sf()
  perimeters <- terra::perim(patchpoly)
  viewshed_areas <- terra::expanse(ptc)
  areas <- terra::expanse(patchpoly)
  total_perimeters <- sum(perimeters)
  total_areas <- sum(areas)
  total_viewshed_areas <- sum(viewshed_areas)
  if (sf::st_crs(m)$units == "ft") {
    num_pt <- round(total_perimeters/3.281)
  } else {
    num_pt <- round(total_perimeters)
  }
  # Number of patches
  Nump <- length(areas)
  # Mean shape index
  MSI <- mean(perimeters/areas)
  # Edge density
  ED <- total_perimeters/total_viewshed_areas
  # Patch size
  PS <- total_areas/Nump
  # Patch density
  PD <- Nump/total_viewshed_areas
  # sample points along the edge of patches
  samples <- sf::st_sample(sf::st_cast(ptc_lines$geometry,
                                       "MULTILINESTRING"),
                           num_pt)
  samples <- sf::st_coordinates(samples)[,-3]
  return(list(Nump, MSI, ED, PS, PD, samples))
}

