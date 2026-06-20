// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

inline bool is_missing(double value) {
  return NumericVector::is_na(value) || !R_finite(value);
}

NumericMatrix raycast_pano(NumericMatrix dsm,
                           double vpt_x, double vpt_y, double vpt_z,
                           int pano_height,
                           int pano_width,
                           double xres,
                           double yres,
                           double ares,
                           double step_size,
                           double max_dist,
                           double sky_value,
                           double sky_threshold,
                           double orientation,
                           bool equirectangular) {
  const int nrow = dsm.nrow();
  const int ncol = dsm.ncol();
  const double PI = 3.141592653589793;
  NumericMatrix pano(pano_height, pano_width);
  pano.fill(sky_value);  // initialize with sky

  for (int y = 0; y < pano_height; ++y) {
    double sin_lat;
    double cos_lat;

    if (equirectangular) {
      // Lu et al. (2020) geo-transform: latitude is linearly sampled over
      // the panorama rows, producing spherical/equirectangular rays.
      const double lat = PI / 2.0 - (y + 0.5) * PI / pano_height;
      sin_lat = std::sin(lat);
      cos_lat = std::cos(lat);
    } else {
      // Cylindrical projection: row position is proportional to tan(latitude).
      const double z_plane = (pano_height / 2.0 - (y + 0.5)) * ares;
      const double norm = std::sqrt(1.0 + z_plane * z_plane);
      sin_lat = z_plane / norm;
      cos_lat = 1.0 / norm;
    }

    for (int x = 0; x < pano_width; ++x) {
      const double lon = x * 2.0 * PI / pano_width - PI + orientation;
      const double sin_lon = std::sin(lon);
      const double cos_lon = std::cos(lon);

      // Unit ray direction
      const double vx = cos_lat * sin_lon;
      const double vy = -cos_lat * cos_lon;
      const double vz = sin_lat;

      bool hit = false;

      for (double step = step_size; step <= max_dist; step += step_size) {
        double px = vpt_x + 0.5 + vx * step / xres;
        double py = vpt_y + 0.5 + vy * step / yres;
        double pz = vpt_z + vz * step;

        int row = static_cast<int>(std::floor(py));
        int col = static_cast<int>(std::floor(px));
        if (row < 0 || col < 0 || row >= nrow || col >= ncol) break;

        double terrain_z = dsm(row, col);
        if (is_missing(terrain_z)) continue;
        // If terrain blocks view
        if (terrain_z - pz > -sky_threshold) {
          pano(y, x) = step;
          hit = true;
          break;
        }
      }

      if (!hit) {
        pano(y, x) = sky_value;  // still sky
      }
    }
  }

  return pano;
}

List raycast_pano_with_semantic(NumericMatrix dsm,
                                NumericMatrix semantic_source,
                                double vpt_x, double vpt_y, double vpt_z,
                                double ground_z,
                                int pano_height,
                                int pano_width,
                                double xres,
                                double yres,
                                double ares,
                                double step_size,
                                double max_dist,
                                double sky_value,
                                double sky_threshold,
                                double orientation,
                                bool equirectangular,
                                bool extrude_semantic) {
  const int nrow = dsm.nrow();
  const int ncol = dsm.ncol();
  const double PI = 3.141592653589793;
  NumericMatrix depth(pano_height, pano_width);
  NumericMatrix semantic_view(pano_height, pano_width);
  depth.fill(sky_value);
  semantic_view.fill(sky_value);

  for (int y = 0; y < pano_height; ++y) {
    double sin_lat;
    double cos_lat;

    if (equirectangular) {
      const double lat = PI / 2.0 - (y + 0.5) * PI / pano_height;
      sin_lat = std::sin(lat);
      cos_lat = std::cos(lat);
    } else {
      const double z_plane = (pano_height / 2.0 - (y + 0.5)) * ares;
      const double norm = std::sqrt(1.0 + z_plane * z_plane);
      sin_lat = z_plane / norm;
      cos_lat = 1.0 / norm;
    }

    for (int x = 0; x < pano_width; ++x) {
      const double lon = x * 2.0 * PI / pano_width - PI + orientation;
      const double vx = cos_lat * std::sin(lon);
      const double vy = -cos_lat * std::cos(lon);
      const double vz = sin_lat;

      for (double step = step_size; step <= max_dist; step += step_size) {
        const double px = vpt_x + 0.5 + vx * step / xres;
        const double py = vpt_y + 0.5 + vy * step / yres;
        const double pz = vpt_z + vz * step;

        const int row = static_cast<int>(std::floor(py));
        const int col = static_cast<int>(std::floor(px));
        if (row < 0 || col < 0 || row >= nrow || col >= ncol) break;

        const double terrain_z = dsm(row, col);
        if (is_missing(terrain_z)) continue;
        const double semantic_z = semantic_source(row, col);

        if (equirectangular && extrude_semantic &&
            !is_missing(semantic_z) && semantic_z > 0.5 &&
            pz >= ground_z - sky_threshold && pz <= terrain_z) {
          depth(y, x) = step;
          semantic_view(y, x) = semantic_z;
          break;
        }

        if (terrain_z - pz > -sky_threshold) {
          depth(y, x) = step;
          semantic_view(y, x) = semantic_z;
          break;
        }
      }
    }
  }

  return List::create(_["depth"] = depth, _["semantic"] = semantic_view);
}

// [[Rcpp::export]]
NumericMatrix dsm_to_pano(NumericMatrix dsm,
                          double vpt_x, double vpt_y, double vpt_z,
                          int pano_height = 256,
                          int pano_width = 512,
                          double xres = 1.0,
                          double yres = 1.0,
                          double ares = 0.0174533,         // ~π/180
                          double step_size = 0.5,
                          double max_dist = 500.0,
                          double sky_value = -1.0,
                          double sky_threshold = 3.0,
                          double orientation = 0.0) {
  return raycast_pano(dsm, vpt_x, vpt_y, vpt_z,
                      pano_height, pano_width, xres, yres,
                      ares, step_size, max_dist,
                      sky_value, sky_threshold, orientation, false);
}

// [[Rcpp::export]]
NumericMatrix dsm_to_equirectangular(NumericMatrix dsm,
                                     double vpt_x, double vpt_y, double vpt_z,
                                     int pano_height = 256,
                                     int pano_width = 512,
                                     double xres = 1.0,
                                     double yres = 1.0,
                                     double step_size = 0.5,
                                     double max_dist = 500.0,
                                     double sky_value = -1.0,
                                     double sky_threshold = 3.0,
                                     double orientation = 0.0) {
  return raycast_pano(dsm, vpt_x, vpt_y, vpt_z,
                      pano_height, pano_width, xres, yres,
                      0.0, step_size, max_dist,
                      sky_value, sky_threshold, orientation, true);
}

// [[Rcpp::export]]
List dsm_semantic_to_pano(NumericMatrix dsm,
                      NumericMatrix semantic,
                      double vpt_x, double vpt_y, double vpt_z,
                      int pano_height = 256,
                      int pano_width = 512,
                      double xres = 1.0,
                      double yres = 1.0,
                      double ares = 0.0174533,
                      double step_size = 0.5,
                      double max_dist = 500.0,
                      double sky_value = -1.0,
                      double sky_threshold = 3.0,
                      double orientation = 0.0) {
  return raycast_pano_with_semantic(dsm, semantic, vpt_x, vpt_y, vpt_z,
                                vpt_z,
                                pano_height, pano_width, xres, yres,
                                ares, step_size,
                                max_dist, sky_value, sky_threshold,
                                orientation, false, false);
}

// [[Rcpp::export]]
List dsm_to_color_equirectangular(NumericMatrix dsm,
                                  NumericMatrix sat_r,
                                  NumericMatrix sat_g,
                                  NumericMatrix sat_b,
                                  double vpt_x, double vpt_y, double vpt_z,
                                  int pano_height = 256,
                                  int pano_width = 512,
                                  double xres = 1.0,
                                  double yres = 1.0,
                                  double step_size = 0.5,
                                  double max_dist = 500.0,
                                  double sky_threshold = 3.0,
                                  double orientation = 0.0,
                                  double sky_r = 135.0,
                                  double sky_g = 206.0,
                                  double sky_b = 235.0) {
  const int nrow = dsm.nrow();
  const int ncol = dsm.ncol();
  const double PI = 3.141592653589793;
  NumericMatrix pano_r(pano_height, pano_width);
  NumericMatrix pano_g(pano_height, pano_width);
  NumericMatrix pano_b(pano_height, pano_width);
  // pre-fill with sky colour
  pano_r.fill(sky_r);
  pano_g.fill(sky_g);
  pano_b.fill(sky_b);

  for (int y = 0; y < pano_height; ++y) {
    const double lat = PI / 2.0 - (y + 0.5) * PI / pano_height;
    const double sin_lat = std::sin(lat);
    const double cos_lat = std::cos(lat);

    for (int x = 0; x < pano_width; ++x) {
      const double lon = x * 2.0 * PI / pano_width - PI + orientation;
      const double vx = cos_lat * std::sin(lon);
      const double vy = -cos_lat * std::cos(lon);
      const double vz = sin_lat;

      for (double step = step_size; step <= max_dist; step += step_size) {
        double px = vpt_x + 0.5 + vx * step / xres;
        double py = vpt_y + 0.5 + vy * step / yres;
        double pz = vpt_z + vz * step;

        int row = static_cast<int>(std::floor(py));
        int col = static_cast<int>(std::floor(px));
        if (row < 0 || col < 0 || row >= nrow || col >= ncol) break;

        double terrain_z = dsm(row, col);
        if (is_missing(terrain_z)) continue;

        if (terrain_z - pz > -sky_threshold) {
          // sample satellite colour at the hit cell
          double r = sat_r(row, col);
          double g = sat_g(row, col);
          double b = sat_b(row, col);
          pano_r(y, x) = is_missing(r) ? sky_r : r;
          pano_g(y, x) = is_missing(g) ? sky_g : g;
          pano_b(y, x) = is_missing(b) ? sky_b : b;
          break;
        }
      }
    }
  }

  return List::create(_["r"] = pano_r, _["g"] = pano_g, _["b"] = pano_b);
}

// [[Rcpp::export]]
List dsm_semantic_to_equirectangular(NumericMatrix dsm,
                                 NumericMatrix semantic,
                                 double vpt_x, double vpt_y, double vpt_z,
                                 double ground_z,
                                 int pano_height = 256,
                                 int pano_width = 512,
                                 double xres = 1.0,
                                 double yres = 1.0,
                                 double step_size = 0.5,
                                 double max_dist = 500.0,
                                 double sky_value = -1.0,
                                 double sky_threshold = 3.0,
                                 double orientation = 0.0,
                                 bool extrude_semantic = false) {
  return raycast_pano_with_semantic(dsm, semantic, vpt_x, vpt_y, vpt_z,
                                ground_z,
                                pano_height, pano_width, xres, yres,
                                0.0, step_size,
                                max_dist, sky_value, sky_threshold,
                                orientation, true, extrude_semantic);
}
