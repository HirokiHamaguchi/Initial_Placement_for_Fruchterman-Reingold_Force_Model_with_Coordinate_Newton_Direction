#pragma once

#include <cmath>
#include <iostream>
#include <utility>

std::pair<float, float> computeDxDy(float gx, float gy, float hxx, float hxy, float hyy,
                                    [[maybe_unused]] float k) {
  // constant step size
  // if (false) {
  //   float stepSizeCoeff = k / std::hypot(gx, gy);
  //   gx *= stepSizeCoeff;
  //   gy *= stepSizeCoeff;
  //   return {-gx, -gy};
  // }

  // positive definite (convex)
  float det = hxx * hyy - hxy * hxy;
  assert(det >= -1e-9);
  float inv_det = 1.0 / det;
  float newton_x = inv_det * (-hyy * gx + hxy * gy);
  float newton_y = inv_det * (-hxx * gy + hxy * gx);
  return {newton_x, newton_y};
}
