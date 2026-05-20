#pragma once

#include "forceModel.hpp"
#include <Eigen/Core>
#include <cassert>
#include <cmath>

/**
 * Fruchterman-Reingold force model
 * Energy: Repulsive = -k^2 * sum(log(d_ij)), Attractive = sum(w_ij * d_ij^3 / (3k))
 */
class FRModel : public ForceModel
{
public:
  FRModel(size_t _n, size_t _m, double _k, std::vector<size_t> &_row,
          std::vector<size_t> &_col, std::vector<double> &_data)
      : k(_k)
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
      double k2 = std::pow(k, 2);
      for (size_t u = 0; u < n; ++u)
      {
        for (size_t v = u + 1; v < n; ++v)
        {
          float d = std::max(1e-10f, (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm());
          score -= k2 * std::log(d);
        }
      }
    }

    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      assert(u < v);
      double w = data[i];
      double d = (position.segment<2>(2 * u) - position.segment<2>(2 * v)).norm();
      score += w * std::pow(d, 3) / (3.0 * k);
    }

    return score;
  }

  double calc_score_and_grad(const Eigen::VectorXf &x,
                             Eigen::VectorXf &grad) const override
  {
    double score = 0.0;
    grad.setZero();

    const double k2 = std::pow(k, 2);
    for (size_t u = 0; u < n; ++u)
    {
      for (size_t v = u + 1; v < n; ++v)
      {
        double dx = x[2 * u] - x[2 * v];
        double dy = x[2 * u + 1] - x[2 * v + 1];
        double d = std::hypot(dx, dy);
        assert(d > 1e-9);
        score -= k2 * std::log(d);

        double g = -k2 / std::pow(d, 2);
        grad[2 * u] += g * dx;
        grad[2 * u + 1] += g * dy;
        grad[2 * v] -= g * dx;
        grad[2 * v + 1] -= g * dy;
      }
    }

    for (size_t i = 0; i < m; ++i)
    {
      size_t u = row[i];
      size_t v = col[i];
      assert(u < v);
      double a = data[i];
      double dx = x[2 * u] - x[2 * v];
      double dy = x[2 * u + 1] - x[2 * v + 1];
      double d = std::hypot(dx, dy);
      score += a * std::pow(d, 3) / (3.0 * k);

      double g = a * d / k;
      grad[2 * u] += g * dx;
      grad[2 * u + 1] += g * dy;
      grad[2 * v] -= g * dx;
      grad[2 * v + 1] -= g * dy;
    }

    return score;
  }

  void calc_grad_hess(float dist, float dx, float dy, float w, float &gx,
                      float &gy, float &hxx, float &hxy, float &hyy) const override
  {
    // Only use attractive force
    float coeff1 = w * dist / k;
    float coeff2 = w / (dist * k);
    gx += coeff1 * dx;
    gy += coeff1 * dy;
    hxx += coeff1 + coeff2 * dx * dx;
    hxy += coeff2 * dx * dy;
    hyy += coeff1 + coeff2 * dy * dy;
  }

  double optimalScaling(Eigen::VectorXf &position) const override
  {
    // Minimize_x x^3 score_a - k^2 n(n-1)/2 \log(x)
    // where score_a = \sum_{i < j} w_{ij} d_{ij}^3 / (3k)
    // Minimize f(x) = x^3 score_a - coeff_r \log(x) : convex
    // f'(x) = 3x^2 score_a - coeff_r / x
    double score_a = calcScore(position, false);
    double coeff_r = std::pow(k, 2) * n * (n - 1) / 2;
    double xStar = std::pow(coeff_r / (3 * score_a), 1.0 / 3);
    position *= xStar;
    return xStar;
  }

private:
  double k; // repulsive constant = 1/sqrt(n)
};
