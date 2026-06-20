#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

struct Vector3 {
  double x, y, z;
};

Vector3 normalize(Vector3 v) {
  const double length = sqrt(pow(v.x, 2) +
                             pow(v.y, 2) +
                             pow(v.z, 2));
  if (length == 0) {
    return {0.0, 0.0, 0.0};
  }
  return {v.x / length, v.y / length, v.z / length};
}

// Compute the outward surface normal directly from the DSM using finite differences.
// Coordinate system: i = column = east (+x), j = row = south (+y), z = up.
// For surface z = f(i,j):  n = normalize( -dz/di, -dz/dj, 1/resolution )
Vector3 cellNormalFromDSM(const Rcpp::NumericMatrix &dsm,
                          int i, int j,
                          double resolution) {
  const int rows = dsm.rows();
  const int cols = dsm.cols();

  // clamp indices so edge cells get a valid 1-sided difference
  int i0 = std::max(i - 1, 0),        i2 = std::min(i + 1, cols - 1);
  int j0 = std::max(j - 1, 0),        j2 = std::min(j + 1, rows - 1);

  double dz_di = (dsm(j, i2) - dsm(j, i0)) / ((i2 - i0) * resolution);
  double dz_dj = (dsm(j2, i) - dsm(j0, i)) / ((j2 - j0) * resolution);

  return normalize({-dz_di, -dz_dj, 1.0});
}

double dotProduct(Vector3 a, Vector3 b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

// double cosAB(int xyp, double zp,
//              int xyt, double zt,
//              double xyn, double zn) {
//   double res;
//   const double pn = sqrt(pow(xyp-xyn, 2) + pow(zp-zn, 2));
//   const double pt = sqrt(pow(xyp-xyt, 2) + pow(zp-zt, 2));
//   const double tn = sqrt(pow(xyn-xyt, 2) + pow(zn-zt, 2));
//   res = (pow(pt,2) + pow(tn,2) - pow(pn,2))/(2*pt*tn);
//   return res;
// }

double PTdistance(int xp, int yp, double zp,
                  int xt, int yt, double zt,
                  double resolution) {
  double res = sqrt(pow((xp-xt)*resolution,2)+
                    pow((yp-yt)*resolution,2)+
                    pow((zp-zt),2));
  return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix VM(const Rcpp::IntegerMatrix &viewshed,
                       const Rcpp::NumericMatrix &dsm,
                       const Rcpp::NumericVector viewpt,
                       const double h,
                       const double resolution) {
  const int rows = dsm.rows();
  const int cols = dsm.cols();
  const double zp = dsm(viewpt[1],viewpt[0]) + h;
  Rcpp::NumericMatrix magnitude(rows, cols);
  Vector3 view = {double(viewpt[0]), double(viewpt[1]), zp};

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      if (i != viewpt[0] || j != viewpt[1]) {
        if (viewshed(j,i) > 0) {
          double zt = dsm(j,i);
          double dis = PTdistance(viewpt[0], viewpt[1], zp,
                                  i, j, zt, resolution);
          // view vector from cell to observer
          Vector3 viewVector = normalize({(view.x - double(i)) * resolution,
                                          (view.y - double(j)) * resolution,
                                          view.z - zt});
          // surface normal from DSM finite differences
          Vector3 normal = cellNormalFromDSM(dsm, i, j, resolution);
          // only front-facing surfaces contribute (dot > 0 means facing observer)
          double angleFactor = std::max(0.0, dotProduct(normal, viewVector));
          if (dis > 0) {
            magnitude(j, i) = angleFactor * resolution*resolution/(dis*dis);
          }
        } else {
          magnitude(j, i) = -9;
        }
      }
    }
  }
  return magnitude;
}
