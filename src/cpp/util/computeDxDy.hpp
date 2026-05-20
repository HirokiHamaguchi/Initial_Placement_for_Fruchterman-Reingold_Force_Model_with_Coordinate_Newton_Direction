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
                                    [[maybe_unused]] float k)
{
  // H = [ hxx hxy ]
  //     [ hxy hyy ]
  // (H + lambda I) dx = -g

  constexpr float eps = 1e-6f;
  float lambda = 1e-6f;
  float a, b, c, det;
  for (int iter = 0; iter < 10; ++iter)
  {
    a = hxx + lambda;
    b = hxy;
    c = hyy + lambda;
    det = a * c - b * b;
    if (a > eps && det > eps)
    {
      break;
    }
    lambda = lambda * 10.0f;
  }

  // Final safeguard. Fallback to gradient descent direction.
  if (!(det > eps))
    return {-gx * 1e-3f, -gy * 1e-3f};

  float inv_det = 1.0f / det;
  float dx = inv_det * (-c * gx + b * gy);
  float dy = inv_det * (-a * gy + b * gx);
  return {dx, dy};
}
