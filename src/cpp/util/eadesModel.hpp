#pragma once

#include "forceModel.hpp"
#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <algorithm>

/**
 * Eades force model
 * Energy: sum(k1 * d_ij * (log(d_ij / a_ij) - 1)) + sum(k2 / (d_ij + eps))
 * Generally more stable than HC model
 */
class EadesModel : public ForceModel
{
public:
  EadesModel(size_t _n, size_t _m, std::vector<size_t> &_row,
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

    // Attractive: sum(k1 * d_ij * (log(d_ij / a_ij) - 1))
    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      double a = data[i]; // target distance
      double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      d = std::max((double)epsilon_r, d);
      score += k1 * d * (std::log(d / a) - 1.0);
    }

    return score;
  }

  double calc_score_and_grad(const Eigen::VectorXf &x,
                             Eigen::VectorXf &grad) const override
  {
    double score = 0.0;
    grad.setZero();

    // Repulsive forces
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

        // gradient: -k2 / (d + eps)^2 * du/dx_i = -k2 / (d + eps)^2 * (x_i - x_j) / d
        double g = -k2 / std::pow(d + epsilon_r, 2);
        grad[2 * u] += g * dx / d;
        grad[2 * u + 1] += g * dy / d;
        grad[2 * v] -= g * dx / d;
        grad[2 * v + 1] -= g * dy / d;
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

      // score: k1 * d * (log(d / a) - 1)
      score += k1 * d * (std::log(d / a) - 1.0);

      // gradient: k1 * log(d / a) * (x_i - x_j) / d
      double g = k1 * std::log(d / a);
      grad[2 * u] += g * dx / d;
      grad[2 * u + 1] += g * dy / d;
      grad[2 * v] -= g * dx / d;
      grad[2 * v + 1] -= g * dy / d;
    }

    return score;
  }

  void calc_grad_hess(float dist, float dx, float dy, float a, float &gx,
                      float &gy, float &hxx, float &hxy, float &hyy) const override
  {
    // Attractive force (discrete phase)
    // gradient: k1 * log(d / a) * u_ij
    // Approximate Hessian: use stabilized version or damped

    float d = std::max((float)epsilon_r, dist);
    float log_ratio = std::log(d / a);
    float coeff_g = k1 * log_ratio / d;

    gx += coeff_g * dx;
    gy += coeff_g * dy;

    // Simplified Hessian: approximate with diagonal dominance for stability
    // Full: H = k1 * [log(d/a) * (I - u*u^T) + 1/d * u*u^T]
    // Stabilized: H ≈ k1 * (log(d/a) / d + 1/d) * I (diagonal)
    // But we need Gauss-Newton-like form: u*u^T
    float coeff_h = k1 / d;
    hxx += coeff_h * dx * dx;
    hxy += coeff_h * dx * dy;
    hyy += coeff_h * dy * dy;
  }

  void optimalScaling(Eigen::VectorXf &position) const override
  {
    // We minimize:
    //   phi(s) = sum k1 * s*d_ij * (log(s*d_ij / a_ij) - 1) + sum k2 / (s*d_ij)
    //          = Cql * s log s + C1 * s + Cm1 / s
    // where Cql = k1 * sum d, C1 = k1 * sum d * (log(d) - log(a) - 1), Cm1 = k2 * sum 1/d

    // Derivatives:
    //   phi'(s) = Cql * (log(s) + 1) + C1 - Cm1 / s^2
    //   phi''(s) = Cql / s + 2 Cm1 / s^3

    double Cql = 0.0;
    double C1 = 0.0;
    double Cm1 = 0.0;

    for (size_t i = 0; i < m; ++i)
    {
      const size_t u = row[i];
      const size_t v = col[i];
      const double a = data[i];
      const double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      Cql += k1 * d;
      C1 += k1 * d * (std::log(d) - std::log(a) - 1.0);
    }
    for (size_t u = 0; u < n; ++u)
    {
      for (size_t v = u + 1; v < n; ++v)
      {
        const double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
        Cm1 += k2 / std::max((double)epsilon_r, d);
      }
    }

    double s = 1.0;
    for (int iter = 0; iter < 20; ++iter)
    {
      const double grad = Cql * (std::log(s) + 1.0) + C1 - Cm1 / (s * s);
      const double hess = Cql / s + 2.0 * Cm1 / (s * s * s);
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
  double k1;        // attractive constant
  double k2;        // repulsive constant
  double epsilon_r; // small value to avoid log(0) and division by zero
};
