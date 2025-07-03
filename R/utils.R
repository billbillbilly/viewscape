#' @import sf
#' @importFrom graphics par
#' @importFrom grDevices rgb
#' @importFrom methods new
#' @importFrom stats sd
#' @importFrom sp SpatialPoints
#' @importFrom torch torch_arange torch_meshgrid nnf_grid_sample torch_stack torch_long
#' @importFrom torch cuda_is_available torch_device torch_zeros_like torch_zeros torch_cat
#' @importFrom terra as.lines

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
  z_viewpoint = terra::extract(dsm,cbind(viewPt[1],viewPt[2]))[1,1]+offset
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

#### geo_process_layer

get_device <- function() {
  if (torch::cuda_is_available()) torch::torch_device("cuda") else torch::torch_device("cpu")
}

# Move tensor to device
to_device <- function(x) {
  x$to(device = get_device())
}



# Normalized grid generator
generate_grid <- function(h, w) {
  x <- torch::torch_linspace(-1.0, 1.0, w)
  y <- torch::torch_linspace(-1.0, 1.0, h)
  xy <- torch::torch_meshgrid(list(x,y), indexing = 'xy')
  grid <- torch::torch_stack(list(xy[[1]], xy[[2]]), dim = -1)
  grid <- grid$unsqueeze(1)$to(dtype = torch::torch_float())
}

# Convert DSM + semantic to voxel
depth2voxel <- function(img_depth, img_sem, gsize, view_row, view_col, view_height) {
  gsize <- torch::torch_tensor(gsize, dtype = torch::torch_int())
  gsize <- as.integer(gsize)
  n <- img_depth$size(1)
  c <- img_depth$size(2)
  h <- img_depth$size(3)
  w <- img_depth$size(4)

  site_z <- img_depth[ , , view_row, view_col] + view_height
  voxel_sitez <- site_z$view(c(n, 1, 1, 1))$expand(c(n, gsize, gsize, gsize))

  # depth voxel
  grid_mask <- generate_grid(gsize, gsize)
  grid_mask <- grid_mask$expand(c(n, gsize, gsize, 2))
  grid_depth <- torch::nnf_grid_sample(img_depth,
                                       grid_mask,
                                       "bilinear",
                                       "zeros",
                                       align_corners=TRUE)

  voxel_depth <- grid_depth$expand(c(n, gsize, gsize, gsize))
  voxel_depth <- voxel_depth - voxel_sitez

  # semantic voxel
  grid_s <- torch::nnf_grid_sample(img_sem$view(c(n, 1, h, w)),
                                   grid_mask,
                                   "nearest",
                                   "zeros",
                                   align_corners=TRUE)
  voxel_s <- grid_s$expand(c(n, gsize, gsize, gsize))

  k <- as.integer(1)
  voxel_s <- grid_s$expand(c(n, gsize/2+k, gsize, gsize))
  gound_s <- torch::torch_zeros(c(n, round(gsize/2-k), gsize, gsize), dtype = torch::torch_float())
  voxel_s <- torch::torch_cat(list(gound_s, voxel_s), 2)

  # occupancy voxel
  voxel_grid <- torch::torch_arange(-gsize/2, gsize/2-1, 1L, dtype = torch::torch_float())
  voxel_grid <- voxel_grid$view(c(1, gsize, 1, 1))$expand(c(n, gsize, gsize, gsize))
  voxel_ocupy <- torch::torch_ge(voxel_depth, voxel_grid)
  voxel_ocupy[ ,gsize-1, , ] <- 0

  # distance voxel
  voxel_dx <- grid_mask[1, , ,1]$view(c(1,1,gsize,gsize))$expand(c(n,gsize,gsize,gsize)) * gsize/2.0
  voxel_dy <- grid_mask[1, , ,2]$view(c(1,1,gsize,gsize))$expand(c(n,gsize,gsize,gsize)) * gsize/2.0
  voxel_dz <- voxel_grid

  voxel_dis <- voxel_dx$mul(voxel_dx) + voxel_dy$mul(voxel_dy) + voxel_dz$mul(voxel_dz)
  voxel_dis <- voxel_dis$add(0.01)   # avoid 1/0 = nan
  voxel_dis <- voxel_dis$mul(voxel_ocupy)
  voxel_dis <- torch::torch_sqrt(voxel_dis) - voxel_ocupy$add(-1.0)$mul(gsize*0.9)

  # Distance computation
  voxel_dx <- grid_mask[1, , , 1]$view(c(1, 1, gsize, gsize))$expand(c(n, gsize, gsize, gsize)) * (gsize / 2)
  voxel_dy <- grid_mask[1, , , 2]$view(c(1, 1, gsize, gsize))$expand(c(n, gsize, gsize, gsize)) * (gsize / 2)
  voxel_dz <- voxel_grid

  # facade:128, tree:64, ground:0, sky:255
  voxel_s = voxel_s$mul(voxel_ocupy) - voxel_ocupy$add(-1.0)$mul(255)

  return(list(voxel_dis = voxel_dis, voxel_s = voxel_s))
}

# Project voxel grid to cylindrical panorama
voxel2pano <- function(voxel_dis, voxel_s, ori, size_pano) {
  PI <- 3.1415926535
  r <- size_pano[1]
  c <- size_pano[2]
  voxel_dis_s <- voxel_dis$size()
  n <- voxel_dis_s[1]
  s <- voxel_dis_s[2]
  t <- voxel_dis_s[3]
  tt <- voxel_dis_s[4]
  k <- as.integer(s/2)

  # rays
  ori <- ori$view(c(n, 1))$expand(c(n, c))
  x <- torch::torch_arange(0, c-1, 1L, dtype = torch::torch_float())$view(c(1, c))$expand(c(n, c))
  y <- torch::torch_arange(0, r-1, 1L, dtype = torch::torch_float())$view(c(1, r))$expand(c(n, r))
  lon <- x * 2 * PI/c + ori - PI
  lat <- PI/2.0 - y * PI/r
  sin_lat <- torch::torch_sin(lat)$view(c(n, 1, r, 1))$expand(c(n, 1, r, c))
  cos_lat <- torch::torch_cos(lat)$view(c(n, 1, r, 1))$expand(c(n, 1, r, c))
  sin_lon <- torch::torch_sin(lon)$view(c(n, 1, 1, c))$expand(c(n, 1, r, c))
  cos_lon <- torch::torch_cos(lon)$view(c(n, 1, 1, c))$expand(c(n, 1, r, c))
  vx <-  cos_lat$mul(sin_lon)
  vy <- -cos_lat$mul(cos_lon)
  vz <-  sin_lat
  vx <- vx$expand(c(n, k, r, c))
  vy <- vy$expand(c(n, k, r, c))
  vz <- vz$expand(c(n, k, r, c))

  voxel_dis <- voxel_dis$contiguous()$view(c(1, n*s*s*s))
  voxel_s <- voxel_s$contiguous()$view(c(1, n*s*s*s))

  # sample voxels along pano-rays
  d_samples <- torch::torch_arange(0, k-1, 1)$view(c(1, k, 1, 1))$expand(c(n, k, r, c))

  samples_x <- vx$mul(d_samples)$add(k)$to(dtype = torch::torch_long())
  samples_y <- vy$mul(d_samples)$add(k)$to(dtype = torch::torch_long())
  samples_z <- vz$mul(d_samples)$add(k)$to(dtype = torch::torch_long())
  #samples_n <- torch::torch_arange(0, n-1, 1)$view(c(n, 1, 1, 1))$expand(c(n, k, r, c))$to(dtype = torch::torch_long())

  samples_x <- torch::torch_clamp(samples_x, 0, s - 1)
  samples_y <- torch::torch_clamp(samples_y, 0, s - 1)
  samples_z <- torch::torch_clamp(samples_z, 0, s - 1)

  #samples_indices <- samples_n$mul(s*s*s)$add(samples_z$mul(s*s))$add(samples_y$mul(s))$add(samples_x)
  samples_indices <- samples_z$mul(s * s)$add(samples_y$mul(s))$add(samples_x)
  # samples_indices <- samples_indices$view(c(1, n*k*r*c))

  samples_indices <- torch::torch_clamp(samples_indices, 0, s * s * s - 1)
  # samples_indices <- samples_indices[1, ]
  samples_indices <- samples_indices$to(dtype = torch::torch_long())

  # get depth pano
  samples_depth <- torch::torch_index_select(voxel_dis, 1, samples_indices)
  samples_depth <- samples_depth$view(c(n, k, r, c))
  min_depth <- torch::torch_min(samples_depth, 1)
  pano_depth <- min_depth[1]
  pano_depth <- pano_depth$view(c(n, 1, r, c))

  # get sem pano
  idx_z <- min_depth[1]
  idx_y <- torch_arange(0, r, 1)$view(c(1, r, 1))$expand(c(n, r, c))
  idx_x <- torch_arange(0, c, 1)$view(c(1, 1, c))$expand(c(n, r, c))
  idx_n <- torch_arange(0, n, 1)$view(c(n, 1, 1))$expand(c(n, r, c))
  idx <- idx_n$mul(k*r*c)$add(idx_z$mul(r*c))$add(idx_y$mul(c))$add(idx_x)$view(c(1, n*r*c))

  samples_s = torch::torch_index_select(voxel_s, 1, samples_indices)
  pano_sem = torch::torch_index_select(samples_s, 1, idx[1, ])$view(c(n,1,r,c))

  return(list(pano_depth = pano_depth, pano_sem = pano_sem))
}

# Full pipeline: DSM + semantic + orientation to pano images
geo_projection <- function(sate_depth, sate_sem, orientation, gsd, pano_size, view_row, view_col, view_height) {
  if (is.null(sate_sem)) {
    sate_sem <- torch::torch_zeros_like(sate_depth)
  }
  # step1: depth to voxel
  gsize <- as.integer(sate_depth$size(4) * gsd)
  voxels <- depth2voxel(sate_depth, sate_sem, gsize, view_row, view_col, view_height)

  pano <- voxel2pano(voxels$voxel_dis, voxels$voxel_s, orientation, pano_size)

  pano$pano_depth <- (pano$pano_depth / 116 - 0.5) / 0.5
  pano$pano_sem <- (pano$pano_sem / 255 - 0.5) / 0.5

  return(pano)
}
