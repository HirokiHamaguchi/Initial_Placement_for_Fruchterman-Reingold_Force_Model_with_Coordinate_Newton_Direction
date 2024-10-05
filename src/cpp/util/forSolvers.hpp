#pragma once

#include <cmath>
#include <iostream>
#include <utility>

std::pair<double, double> computeEigenvalues(double h_xx, double h_xy, double h_yy) {
  // The characteristic equation is of the form: λ^2 - trace * λ + determinant = 0
  // Solving this quadratic equation: λ = (trace ± sqrt(trace^2 - 4 * determinant)) / 2

  double trace = h_xx + h_yy;
  double determinant = h_xx * h_yy - h_xy * h_xy;
  double discriminant = trace * trace - 4 * determinant;
  assert(discriminant >= -1e-9);  // The symmetric matrix should have real eigenvalues

  double lambda1 = (trace + std::sqrt(discriminant + 1e-9)) / 2.0;
  double lambda2 = (trace - std::sqrt(discriminant + 1e-9)) / 2.0;

  return {lambda1, lambda2};
}

std::pair<double, double> computeDxDy(double gx, double gy, double hxx, double hxy,
                                      double hyy, [[maybe_unused]] double k) {
  auto [lambda1, lambda2] = computeEigenvalues(hxx, hxy, hyy);
  if (lambda1 > 1e-2 && lambda2 > 1e-2) {
    // positive definite (convex quadratic)
    double det = hxx * hyy - hxy * hxy;
    assert(det >= -1e-9);
    double inv_det = 1.0 / det;
    double newton_x = inv_det * (-hyy * gx + hxy * gy);
    double newton_y = inv_det * (-hxx * gy + hxy * gx);
    return {newton_x, newton_y};
  } else {
    assert(false);
  }
}
