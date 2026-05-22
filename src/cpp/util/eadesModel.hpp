#pragma once

#include "forceModel.hpp"
#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "dbg.h"

/**
 * Eades force model
 * Energy: sum(k1 * d_ij * (log(d_ij / a_ij) - 1)) + sum(k2 / (d_ij))
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

  double calcScore(const Eigen::VectorXd &position,
                   bool includeRepulsive = true) const override
  {
    double score = 0.0;

    if (includeRepulsive)
    {
      // Repulsive: sum(k2 / (d_ij))
      for (size_t u = 0; u < n; ++u)
      {
        for (size_t v = u + 1; v < n; ++v)
        {
          double d = std::max(epsilon_r,
                              (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());
          score += k2 / d;
        }
      }
    }

    // Attractive: sum(k1 * d_ij * (log(d_ij / a_ij) - 1))
    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      double a = std::max(epsilon_r, data[i]); // target distance
      double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      d = std::max(epsilon_r, d);
      score += k1 * d * (std::log(d / a) - 1.0);
    }

    return score;
  }

  double calc_score_and_grad(const Eigen::VectorXd &x,
                             Eigen::VectorXd &grad) const override
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
        double d = std::max(epsilon_r, std::hypot(dx, dy));
        score += k2 / d;

        // gradient: -k2 / d^2 * (x_i - x_j) / d = -k2 * (x_i - x_j) / d^3
        double g = -k2 / std::pow(d, 2);
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
      double dx = x[2 * u] - x[2 * v];
      double dy = x[2 * u + 1] - x[2 * v + 1];
      double a = std::max(epsilon_r, data[i]); // target distance
      double d = std::max(epsilon_r, std::hypot(dx, dy));
      double log_da = std::log(d / a);

      // score: k1 * d * (log(d / a) - 1)
      score += k1 * d * (log_da - 1.0);

      // gradient: k1 * log(d / a) * (x_i - x_j) / d
      double g = k1 * log_da;
      grad[2 * u] += g * dx / d;
      grad[2 * u + 1] += g * dy / d;
      grad[2 * v] -= g * dx / d;
      grad[2 * v + 1] -= g * dy / d;
    }

    return score;
  }

  void calc_grad_hess(double dist, double dx, double dy, double a, double &gx,
                      double &gy, double &hxx, double &hxy, double &hyy) const override
  {
    // Attractive force (discrete phase)
    // gradient: k1 * log(d / a) * (x_i - x_j) / d
    double d = std::max(epsilon_r, dist);
    double log_ratio = std::log(d / a);
    double coeff_g = k1 * log_ratio / d;
    gx += coeff_g * dx;
    gy += coeff_g * dy;

    // H = k1/d * (log(d/a) * I + (1 - log(d/a)) * (x_i - x_j)(x_i - x_j)^T / d^2)
    double coeff_h = k1 / d;
    hxx += coeff_h * (log_ratio + (1.0 - log_ratio) * (dx * dx) / (d * d));
    hxy += coeff_h * ((1.0 - log_ratio) * (dx * dy) / (d * d));
    hyy += coeff_h * (log_ratio + (1.0 - log_ratio) * (dy * dy) / (d * d));
  }

  double optimalScaling(Eigen::VectorXd &position) const override
  {
    // We minimize:
    //   phi(s) = sum k1 * s*d_ij * (log(s*d_ij / a_ij) - 1) + sum k2 / (s*d_ij)
    //          = Cql * s log s + C1 * s + Cm1 / s
    //
    // Let s = exp(t), then:
    //
    //   psi(t) = phi(exp(t))
    //          = Cql * exp(t) * t + C1 * exp(t) + Cm1 * exp(-t)
    //
    // Derivatives:
    //
    //   psi'(t)  = Cql * exp(t) * (t + 1)
    //             + C1 * exp(t)
    //             - Cm1 * exp(-t)
    //
    //   psi''(t) = Cql * exp(t) * (t + 2)
    //             + C1 * exp(t)
    //             + Cm1 * exp(-t)

    double Cql = 0.0;
    double C1 = 0.0;
    double Cm1 = 0.0;

    for (size_t i = 0; i < m; ++i)
    {
      const size_t u = row[i];
      const size_t v = col[i];
      const double a = data[i];
      const double d = std::max(
          epsilon_r,
          (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());

      Cql += d;
      C1 += d * (std::log(d / a) - 1.0);
    }

    Cql *= k1;
    C1 *= k1;

    for (size_t u = 0; u < n; ++u)
    {
      for (size_t v = u + 1; v < n; ++v)
      {
        const double d = std::max(
            epsilon_r,
            (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());

        Cm1 += k2 / d;
      }
    }

    // Optimize in log-space: s = exp(t)
    double t = std::log(0.1);

    for (int iter = 0; iter < 20; ++iter)
    {
      const double et = std::exp(t);
      const double emt = 1.0 / et;

      const double grad =
          Cql * et * (t + 1.0) + C1 * et - Cm1 * emt;

      const double hess =
          Cql * et * (t + 2.0) + C1 * et + Cm1 * emt;

      const double step = grad / hess;
      const double t_new = t - step;

      if (std::abs(t_new - t) < 1e-8)
      {
        t = t_new;
        break;
      }

      t = t_new;
    }

    const double s = std::exp(t);

    position *= s;
    return s;
  }

private:
  double k1;        // attractive constant
  double k2;        // repulsive constant
  double epsilon_r; // small value to avoid log(0) and division by zero
};
