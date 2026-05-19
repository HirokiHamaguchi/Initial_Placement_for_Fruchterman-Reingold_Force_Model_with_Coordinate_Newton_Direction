#pragma once

#include "forceModel.hpp"
#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <algorithm>

/**
 * Hooke-Coulomb force model
 * Energy (attractive): sum(k1/2 * (d_ij - a_ij)^2) where a_ij = data[i]
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
      double a = data[i]; // target distance
      double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      double diff = d - a;
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
    // phi(s) = sum(k1/2 * (s*d_ij - a_ij)^2) + sum(k2 / (s*d_ij + eps))
    // phi'(s) = k1 * sum(d_ij * (s*d_ij - a_ij)) - k2 * sum(d_ij / (s*d_ij + eps)^2)
    // Solve phi'(s) = 0 numerically using Newton's method

    double s = 1.0;
    double score_r = 0.0;

    // Compute repulsive energy contribution
    for (size_t u = 0; u < n; ++u)
    {
      for (size_t v = u + 1; v < n; ++v)
      {
        float d = std::max((float)epsilon_r,
                           (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());
        score_r += k2 / (d + epsilon_r);
      }
    }

    // Newton's method: s_{n+1} = s_n - f(s_n) / f'(s_n)
    for (int iter = 0; iter < 10; ++iter)
    {
      double phi = 0.0, phi_prime = 0.0;

      // Attractive contribution
      for (size_t i = 0; i < m; ++i)
      {
        size_t u = row[i];
        size_t v = col[i];
        double a = data[i];
        double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
        d = std::max((double)epsilon_r, d);
        double sd = s * d;
        double diff = sd - a;
        phi += k1 * 0.5 * diff * diff;
        phi_prime += k1 * d * diff;
      }

      // Repulsive contribution
      for (size_t u = 0; u < n; ++u)
      {
        for (size_t v = u + 1; v < n; ++v)
        {
          double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
          d = std::max(epsilon_r, d);
          double sd = s * d;

          phi += k2 / (sd + epsilon_r);
          phi_prime -= k2 * d / std::pow(sd + epsilon_r, 2);
        }
      }

      if (std::abs(phi_prime) < 1e-12)
        break;

      double s_new = s - phi / phi_prime;
      s_new = std::max(0.1, std::min(10.0, s_new)); // bound to avoid numerical issues

      if (std::abs(s_new - s) < 1e-6)
        break;

      s = s_new;
    }

    position *= s;
  }

private:
  double k1;        // attractive spring constant
  double k2;        // repulsive constant
  double epsilon_r; // small value to avoid division by zero
};
