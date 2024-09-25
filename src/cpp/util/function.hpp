#pragma once

#include <Eigen/Core>

#include "problem.hpp"

struct Function {
  const Problem& problem;
  Function(const Problem& problem) : problem(problem) {}

  float operator()(const Eigen::VectorXf& x, Eigen::VectorXf& grad) {
    double score = 0.0;
    grad.setZero();

    double k2 = std::pow(problem.k, 2);
    for (size_t u = 0; u < problem.n; ++u) {
      for (size_t v = u + 1; v < problem.n; ++v) {
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

    for (size_t i = 0; i < problem.m; ++i) {
      size_t u = problem.row[i];
      size_t v = problem.col[i];
      assert(u < v);
      double a = problem.data[i];
      double dx = x[2 * u] - x[2 * v];
      double dy = x[2 * u + 1] - x[2 * v + 1];
      double d = std::hypot(dx, dy);
      score += a * std::pow(d, 3) / (3.0 * problem.k);

      double g = a * d / problem.k;
      grad[2 * u] += g * dx;
      grad[2 * u + 1] += g * dy;
      grad[2 * v] -= g * dx;
      grad[2 * v + 1] -= g * dy;
    }

    return score;
  }
};