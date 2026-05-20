#pragma once

#include "forceModel.hpp"
#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <algorithm>

/**
 * Hooke-Coulomb force model
 * Energy: sum(k1/2 * (d_ij - a_ij)^2) + sum(k2 / (d_ij + epsilon_r))
 * Uses Gauss-Newton Hessian for numerical stability
 */
class HCModel : public ForceModel
{
public:
  HCModel(size_t _n, size_t _m, std::vector<size_t> &_row,
          std::vector<size_t> &_col, std::vector<double> &_data)
      : k1(1.0), k2(1.0), epsilon_r(1e-10)
  {
    n = _n;
    m = _m;
    row = _row;
    col = _col;
    data = _data;
    assert(row.size() == m);
    assert(col.size() == m);
    assert(data.size() == m);
    makeAdj();
  }

  double calcScore(const Eigen::VectorXf &position,
                   bool includeRepulsive = true) const override
  {
    double score = 0.0;

    if (includeRepulsive)
    {
      // Repulsive: sum(k2 / (d_ij + epsilon_r))
      for (size_t u = 0; u < n; ++u)
      {
        for (size_t v = u + 1; v < n; ++v)
        {
          float d = std::max((float)epsilon_r,
                             (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());
          score += k2 / (d + epsilon_r);
        }
      }
    }

    // Attractive: sum(k1/2 * (d_ij - a_ij)^2)
    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      assert(u < v);
      double w = data[i]; // target distance
      double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      double diff = d - w;
      score += k1 * 0.5 * diff * diff;
    }

    return score;
  }

  double calc_score_and_grad(const Eigen::VectorXf &x,
                             Eigen::VectorXf &grad) const override
  {
    double score = 0.0;
    grad.setZero();

    // Repulsive forces
    if (true) // always include for LBFGS
    {
      for (size_t u = 0; u < n; ++u)
      {
        for (size_t v = u + 1; v < n; ++v)
        {
          double dx = x[2 * u] - x[2 * v];
          double dy = x[2 * u + 1] - x[2 * v + 1];
          double d = std::hypot(dx, dy);
          d = std::max((double)epsilon_r, d);

          // score: k2 / (d + eps)
          score += k2 / (d + epsilon_r);

          // gradient of k2 / (d + eps) w.r.t. d is -k2 / (d + eps)^2
          double g = -k2 / std::pow(d + epsilon_r, 2);
          if (d > epsilon_r)
          {
            grad[2 * u] += g * dx / d;
            grad[2 * u + 1] += g * dy / d;
            grad[2 * v] -= g * dx / d;
            grad[2 * v + 1] -= g * dy / d;
          }
        }
      }
    }

    // Attractive forces
    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      double a = data[i]; // target distance
      double dx = x[2 * u] - x[2 * v];
      double dy = x[2 * u + 1] - x[2 * v + 1];
      double d = std::hypot(dx, dy);
      d = std::max((double)epsilon_r, d);

      // score: k1/2 * (d - a)^2
      double diff = d - a;
      score += k1 * 0.5 * diff * diff;

      // gradient: k1 * (d - a) * u_ij where u_ij = (x_i - x_j) / d
      double g = k1 * diff;
      if (d > epsilon_r)
      {
        grad[2 * u] += g * dx / d;
        grad[2 * u + 1] += g * dy / d;
        grad[2 * v] -= g * dx / d;
        grad[2 * v + 1] -= g * dy / d;
      }
    }

    return score;
  }

  void calc_grad_hess(float dist, float dx, float dy, float a, float &gx,
                      float &gy, float &hxx, float &hxy, float &hyy) const override
  {
    // Attractive force only (for discrete phase)
    // gradient: k1 * (dist - a) * u_ij
    // Gauss-Newton Hessian: k1 * u_ij * u_ij^T

    float d = std::max((float)epsilon_r, dist);
    float diff = d - a;
    float coeff_g = k1 * diff / d; // gradient coefficient
    float coeff_h = k1 / d;        // Hessian coefficient (Gauss-Newton)

    gx += coeff_g * dx;
    gy += coeff_g * dy;
    hxx += coeff_h * dx * dx;
    hxy += coeff_h * dx * dy;
    hyy += coeff_h * dy * dy;
  }

  void optimalScaling(Eigen::VectorXf &position) const override
  {
    //
    // We minimize:
    //   phi(s) = sum_{E}  k1/2 * (s d_ij - a_ij)^2 + sum_{u<v}  k2 / (s d_ij)
    //          = C2 * s^2 + C1 * s + Cm1 / s + const
    // with C2 = (k1/2) * sum d^2, C1 = -k1 * sum d a, Cm1 = k2 * sum 1/d
    //
    // The derivatives are:
    //   phi'(s) = 2 C2 s + C1 - Cm1 / s^2
    //   phi''(s) = 2 C2 + 2 Cm1 / s^3

    double C2 = 0.0;
    double C1 = 0.0;
    double Cm1 = 0.0;

    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      const double a = data[i];
      const double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      C2 += 0.5 * k1 * d * d;
      C1 += -k1 * d * a;
    }
    for (size_t u = 0; u < n; ++u)
    {
      for (size_t v = u + 1; v < n; ++v)
      {
        const double d =
            (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
        Cm1 += k2 / std::max((double)epsilon_r, d);
      }
    }

    // Newton optimization
    double s = 1.0;

    for (int iter = 0; iter < 20; ++iter)
    {
      const double grad = 2.0 * C2 * s + C1 - Cm1 / (s * s);
      const double hess = 2.0 * C2 + 2.0 * Cm1 / (s * s * s);
      const double step = grad / hess;
      double s_new = s - step;
      s_new = std::max(1e-8, s_new);
      if (std::abs(s_new - s) < 1e-8)
      {
        s = s_new;
        break;
      }
      s = s_new;
    }

    position *= s;
  }

private:
  double k1;        // attractive spring constant
  double k2;        // repulsive constant
  double epsilon_r; // small value to avoid division by zero
};
