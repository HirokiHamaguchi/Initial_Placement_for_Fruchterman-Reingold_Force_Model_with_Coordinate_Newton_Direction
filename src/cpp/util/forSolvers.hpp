#pragma once

#include <cmath>
#include <iostream>
#include <utility>

std::pair<double, double> computeDxDy(double gx, double gy, double hxx, double hxy,
                                      double hyy, [[maybe_unused]] double k) {
  // positive definite (convex)
  double det = hxx * hyy - hxy * hxy;
  assert(det >= -1e-9);
  double inv_det = 1.0 / det;
  double newton_x = inv_det * (-hyy * gx + hxy * gy);
  double newton_y = inv_det * (-hxx * gy + hxy * gx);
  return {newton_x, newton_y};
}
