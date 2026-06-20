#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerMatrix sectorMask(
    const Rcpp::IntegerMatrix &viewshed,
    const Rcpp::NumericVector viewpt,
    const Rcpp::NumericVector fov) {
  const int rows = viewshed.rows();
  const int cols = viewshed.cols();
  Rcpp::IntegerMatrix visible(rows, cols);
  const double rad2deg = 180.0 / M_PI;

  const double xp = viewpt[0];  // viewpoint column (east)
  const double yp = viewpt[1];  // viewpoint row

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      if (viewshed(j, i) == 1) {
        double dx = double(i) - xp;         // east component
        double dy = yp - double(j);         // north component (flip row axis)
        // atan2(dx, dy): clockwise from north, range (-180, 180]
        double angle = std::atan2(dx, dy) * rad2deg;

        bool in_sector;
        if (fov[0] <= fov[1]) {
          // normal sector: e.g. c(-60, 60) for a north-facing arc
          in_sector = (angle >= fov[0] && angle <= fov[1]);
        } else {
          // wrap-around crossing ±180 (south): e.g. c(135, -135)
          in_sector = (angle >= fov[0] || angle <= fov[1]);
        }
        if (in_sector) {
          visible(j, i) = viewshed(j, i);
        }
      }
    }
  }
  return visible;
}
