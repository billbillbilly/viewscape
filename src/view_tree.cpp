#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>

using namespace Rcpp;

// Curvature and refraction correction (identical formula used in visibleLabel.cpp)
static inline double curvRefCorr(double distance, double refraction_factor) {
  return (distance * distance) / 12740000.0 * refraction_factor;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix view_tree(
    const Rcpp::NumericVector& viewpoint,
    const Rcpp::NumericMatrix& dsm,
    const double h,
    const int max_dis,
    const double refraction_factor) {

  const int rows = dsm.rows();
  const int cols = dsm.cols();

  // Observer position (col = x, row = y), 0-indexed
  const int col_O = static_cast<int>(viewpoint[0]);
  const int row_O = static_cast<int>(viewpoint[1]);
  const double z_O = viewpoint[2];

  // Buffer radius R (in cells) for the view-block relationship.
  // Wang et al. (2022) recommend R in [0.5*L, L] where L = cell width.
  // With cell units L = 1, R = 0.5 is the tightest (most accurate) setting.
  const double R = 0.5;

  Rcpp::IntegerMatrix VS(rows, cols);           // output: 1 = visible, 0 = not
  std::vector<bool>   VA(rows * cols, false);   // visited flag

  auto in_bounds = [&](int r, int c) -> bool {
    return r >= 0 && r < rows && c >= 0 && c < cols;
  };

  // Returns true if C (col_C, row_C) generates view-block between O and B (col_B, row_B).
  // Implements Wang et al. (2022) equations 7–9:
  //   Eq 7: perpendicular distance from C to line OB <= R
  //   Eq 8: dot(OC, OB) > 0  (C lies on the O-side of B)
  //   Eq 9: dot(BC, BO) > 0  (C lies on the B-side of O)
  auto generates_viewblock = [&](int col_C, int row_C, int col_B, int row_B) -> bool {
    const double dX = col_B - col_O;
    const double dY = row_B - row_O;
    const double dAB_sq = dX * dX + dY * dY;
    if (dAB_sq == 0.0) return false;
    // Eq 7
    const double cross = std::abs((col_C - col_O) * dY - (row_C - row_O) * dX);
    if (cross > R * std::sqrt(dAB_sq)) return false;
    // Eq 8
    if ((col_C - col_O) * dX + (row_C - row_O) * dY <= 0.0) return false;
    // Eq 9
    if ((col_C - col_B) * (-dX) + (row_C - row_B) * (-dY) <= 0.0) return false;
    return true;
  };

  // Mark observer as visible and visited
  VS(row_O, col_O) = 1;
  VA[row_O * cols + col_O] = true;

  // 8-directional neighbor offsets: (drow, dcol)
  constexpr int nbr_dr[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
  constexpr int nbr_dc[8] = {-1,  0,  1, -1, 1, -1, 0, 1};

  // Iterative DFS — stack stores (col, row, mvs)
  struct Item { int col, row; double mvs; };
  std::vector<Item> stack;
  stack.reserve(static_cast<size_t>(rows * cols));

  // Seed with the observer's 8 immediate neighbors
  for (int i = 0; i < 8; ++i) {
    int nr = row_O + nbr_dr[i];
    int nc = col_O + nbr_dc[i];
    if (in_bounds(nr, nc)) {
      stack.push_back({nc, nr, -std::numeric_limits<double>::infinity()});
    }
  }

  while (!stack.empty()) {
    Item cur = stack.back();
    stack.pop_back();

    if (VA[cur.row * cols + cur.col]) continue;

    const double dc = static_cast<double>(cur.col - col_O);
    const double dr = static_cast<double>(cur.row - row_O);
    const double dist_cells = std::sqrt(dc * dc + dr * dr);

    if (dist_cells > static_cast<double>(max_dis)) {
      VA[cur.row * cols + cur.col] = true;
      continue;
    }

    // Target elevation with height offset and curvature/refraction correction.
    // The distance argument follows the same convention as visibleLabel.cpp
    // (dist_cells * h approximates metres when h ~ cell resolution in metres).
    const double dist_m   = dist_cells * h;
    const double z_T      = dsm(cur.row, cur.col) + h - curvRefCorr(dist_m, refraction_factor);

    // Actual slope of the line-of-sight from O to T (Wang et al. 2022, Sec 2.3)
    const double AS       = (z_T - z_O) / dist_cells;

    double mvs_next;
    if (AS > cur.mvs) {
      VS(cur.row, cur.col) = 1;
      mvs_next = AS;
    } else {
      mvs_next = cur.mvs;   // horizon MVS propagates unchanged to children
    }
    VA[cur.row * cols + cur.col] = true;

    // Push child nodes: neighbors of T that T view-blocks (from O)
    for (int i = 0; i < 8; ++i) {
      const int nr = cur.row + nbr_dr[i];
      const int nc = cur.col + nbr_dc[i];
      if (in_bounds(nr, nc) && !VA[nr * cols + nc]) {
        if (generates_viewblock(cur.col, cur.row, nc, nr)) {
          stack.push_back({nc, nr, mvs_next});
        }
      }
    }
  }

  return VS;
}
