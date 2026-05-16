#pragma once

#include <Eigen/Core>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "../util/problem.hpp"
#include "../util/timer.hpp"

static double calcFRScoreMatrix(const Problem &problem, const Eigen::MatrixXf &pos)
{
  double score = 0.0;
  const double k2 = problem.k * problem.k;

  for (size_t u = 0; u < problem.n; ++u)
  {
    for (size_t v = u + 1; v < problem.n; ++v)
    {
      double dx = static_cast<double>(pos(0, u) - pos(0, v));
      double dy = static_cast<double>(pos(1, u) - pos(1, v));
      double d = std::hypot(dx, dy);
      d = std::max(d, 1e-9);
      score -= k2 * std::log(d);
    }
  }

  for (size_t i = 0; i < problem.m; ++i)
  {
    size_t u = problem.row[i];
    size_t v = problem.col[i];

    double a = problem.data[i];
    double dx = static_cast<double>(pos(0, u) - pos(0, v));
    double dy = static_cast<double>(pos(1, u) - pos(1, v));
    double d = std::hypot(dx, dy);
    score += a * std::pow(d, 3) / (3.0 * problem.k);
  }

  return score;
}

void solve_FR(const Problem &problem, std::vector<Eigen::VectorXf> &positions,
              std::vector<std::pair<double, double>> &hist, Timer &timer,
              const int MAX_ITER)
{
  timer.start();

  size_t nnodes = problem.n;
  size_t dim = 2;
  int iterations = MAX_ITER;
  double threshold = 1e-5; // changed from 1e-3
  assert(!positions.empty());

  Eigen::VectorXf posVec = positions.back();
  problem.optimalScaling(posVec);
  Eigen::MatrixXf pos =
      Eigen::Map<const Eigen::MatrixXf>(posVec.data(), dim, nnodes);

  const double k = problem.k;
  const float scaling = pos.array().abs().maxCoeff();
  double t = 0.1 * scaling;
  double dt = t / (iterations + 1);

  Eigen::MatrixXf displacement = Eigen::MatrixXf::Zero(dim, nnodes);
  for (int iteration = 0; iteration < iterations; ++iteration)
  {
    displacement.setZero();

    for (size_t i = 0; i < nnodes; ++i)
    {
      Eigen::MatrixXf delta = pos.col(i).replicate(1, nnodes) - pos;
      Eigen::RowVectorXf distance2 = delta.colwise().squaredNorm();
      distance2 =
          distance2.unaryExpr([scaling](float d)
                              { return std::max(d, 0.01f * 0.01f * scaling * scaling); });
      Eigen::RowVectorXf distance = distance2.array().sqrt();
      displacement.col(i) += (delta.array().rowwise() * (k * k / distance2.array()))
                                 .rowwise()
                                 .sum()
                                 .matrix();
      for (auto &[j, w] : problem.adj[i])
        displacement.col(i) -= delta.col(j) * (w * distance[j] / k);
    }

    Eigen::RowVectorXf length = displacement.colwise().norm();
    length = length.unaryExpr([scaling](float l)
                              { return std::max(l, 0.1f * scaling); });
    Eigen::MatrixXf direction = displacement.array().rowwise() / length.array();

    if (iteration == 0)
    {
      const double current_score = calcFRScoreMatrix(problem, pos);
      while (true)
      {
        Eigen::MatrixXf trial_pos = pos + direction * static_cast<float>(t);
        double trial_score = calcFRScoreMatrix(problem, trial_pos);
        if (trial_score <= current_score || t <= 1e-2 * scaling)
          break;
        t *= 0.5;
      }
      dt = t / (iterations + 1);
    }

    Eigen::MatrixXf delta_pos = direction * static_cast<float>(t);

    pos += delta_pos;
    t -= dt;

    timer.stop();
    positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), pos.size()));
    hist.emplace_back(problem.calcScore(positions.back()), timer.sec());
    timer.start();

    if (delta_pos.norm() / nnodes < threshold * scaling)
      break;
  }

  timer.stop();
}
