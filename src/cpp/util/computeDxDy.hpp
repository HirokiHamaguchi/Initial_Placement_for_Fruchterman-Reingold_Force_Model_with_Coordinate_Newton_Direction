#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

std::pair<float, float> computeDxDy(float gx,
                                    float gy,
                                    float hxx,
                                    float hxy,
                                    float hyy,
                                    float grid_size)
{
  // H = [ hxx hxy ]
  //     [ hxy hyy ]
  float a = hxx, b = hxy, c = hyy, det = a * c - b * b;

  if (det < 1e-9f)
  {
    // If Hessian is nearly singular, use gradient descent with a small step size
    float norm_g = std::sqrt(gx * gx + gy * gy);
    if (norm_g < 1e-9f)
      return {0.0f, 0.0f};
    float step_size = 2.0f * grid_size / norm_g;
    return {-step_size * gx, -step_size * gy};
  }

  float inv_det = 1.0f / det;
  float dx = inv_det * (-c * gx + b * gy);
  float dy = inv_det * (-a * gy + b * gx);
  return {dx, dy};
}
