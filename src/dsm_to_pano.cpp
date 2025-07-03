// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dsm_to_pano(NumericMatrix dsm,
                          double vpt_x, double vpt_y, double vpt_z,
                          int pano_height = 256,
                          int pano_width = 512,
                          double ares = 0.0174533,         // ~π/180
                          double step_size = 0.5,
                          double max_dist = 500.0,
                          double sky_value = -1.0,
                          double sky_threshold = 3.0,
                          double orientation = 0.0) {

  int nrow = dsm.nrow();
  int ncol = dsm.ncol();

  NumericMatrix pano(pano_height, pano_width);
  pano.fill(sky_value);  // initialize with sky

  const double PI = 3.141592653589793;

  for (int y = 0; y < pano_height; ++y) {
    double lat = PI / 2.0 - y * PI / pano_height;
    double sin_lat = std::sin(lat);
    double cos_lat = std::cos(lat);

    for (int x = 0; x < pano_width; ++x) {
      double lon = x * 2.0 * PI / pano_width - PI + orientation;
      double sin_lon = std::sin(lon);
      double cos_lon = std::cos(lon);

      // Unit ray direction
      double vx = cos_lat * sin_lon;
      double vy = -cos_lat * cos_lon;
      double vz = sin_lat;

      bool hit = false;

      for (double step = step_size; step <= max_dist; step += step_size) {
        double px = vpt_x + 0.5 + vx * step;
        double py = vpt_y + 0.5 + vy * step;
        double pz = vpt_z + vz * step;

        int row = static_cast<int>(std::round(py));
        int col = static_cast<int>(std::round(px));
        if (row < 0 || col < 0 || row >= nrow || col >= ncol) break;

        double terrain_z = dsm(row, col);
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
