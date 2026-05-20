#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include "dbg.h"

std::pair<double, double> computeDxDy(double gx,
                                      double gy,
                                      double hxx,
                                      double hxy,
                                      double hyy,
                                      double grid_size)
{
  // H = [ hxx hxy ]
  //     [ hxy hyy ]
  double a = hxx, b = hxy, c = hyy, det = a * c - b * b;
  if (det < 1e-9)
  {
    // If Hessian is nearly singular, use gradient descent with a small step size
    double norm_g = std::sqrt(gx * gx + gy * gy);
    if (norm_g < 1e-9)
      return {0.0, 0.0};
    double step_size = 2.0 * grid_size / norm_g;
    return {-step_size * gx, -step_size * gy};
  }
  double inv_det = 1.0 / det;
  double dx = inv_det * (-c * gx + b * gy);
  double dy = inv_det * (-a * gy + b * gx);
  double norm_d = std::hypot(dx, dy);
  if (norm_d > 3.0 * grid_size)
  {
    double scale = (3.0 * grid_size) / norm_d;
    dx *= scale;
    dy *= scale;
  }
  return {dx, dy};
}
