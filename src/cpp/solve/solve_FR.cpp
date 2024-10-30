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

#include "../util/dbg.h"
#include "../util/problem.hpp"

std::vector<Eigen::VectorXf> solve_FR(const Problem& problem,
                                      std::vector<Eigen::VectorXf>& positions,
                                      std::vector<std::pair<double, double>>& hist,
                                      Timer& timer, const int MAX_ITER) {
  timer.start();

  size_t nnodes = problem.n;
  size_t dim = 2;
  int iterations = MAX_ITER;
  double threshold = 1e-4;  // changed from 1e-3
  assert(!positions.empty());

  Eigen::MatrixXf pos =
      Eigen::Map<const Eigen::MatrixXf>(positions.back().data(), dim, nnodes);

  // rescale layout to fit in a unit square (from networkX rescale_layout)
  pos.colwise() -= pos.rowwise().mean();
  float lim = pos.array().abs().maxCoeff();
  if (lim > 0) pos *= 1.0 / lim;

  const double k = problem.k;  // std::sqrt(1.0 / nnodes);
  double t = 0.1;
  double dt = t / (iterations + 1);

  Eigen::MatrixXf displacement = Eigen::MatrixXf::Zero(dim, nnodes);
  for (int iteration = 0; iteration < iterations; ++iteration) {
    displacement.setZero();

    for (size_t i = 0; i < nnodes; ++i) {
      Eigen::MatrixXf delta = pos.col(i).replicate(1, nnodes) - pos;
      Eigen::RowVectorXf distance2 = delta.colwise().squaredNorm();
      distance2 =
          distance2.unaryExpr([](float d) { return std::max(d, 0.01f * 0.01f); });
      Eigen::RowVectorXf distance = distance2.array().sqrt();
      displacement.col(i) += (delta.array().rowwise() * (k * k / distance2.array()))
                                 .rowwise()
                                 .sum()
                                 .matrix();
      for (auto& [j, w] : problem.adj[i])
        displacement.col(i) -= delta.col(j) * (w * distance[j] / k);
    }

    Eigen::RowVectorXf length = displacement.colwise().norm();
    length = length.unaryExpr([](float l) { return std::max(l, 0.1f); });
    Eigen::MatrixXf delta_pos = displacement.array().rowwise() * (t / length.array());

    pos += delta_pos;
    t -= dt;

    timer.stop();
    positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), pos.size()));
    hist.emplace_back(problem.calcScore(positions.back()), timer.sec());
    timer.start();

    if (delta_pos.norm() / nnodes < threshold) break;
  }

  timer.stop();

  return positions;
}
