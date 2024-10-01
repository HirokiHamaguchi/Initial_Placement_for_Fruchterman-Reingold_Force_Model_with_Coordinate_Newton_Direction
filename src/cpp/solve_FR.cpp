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

#include "util/dbg.h"
#include "util/problem.hpp"

std::vector<Eigen::VectorXf> solve_FR(const Problem& problem,
                                      const Eigen::VectorXf& _pos, int iterations,
                                      double threshold) {
  std::vector<Eigen::VectorXf> positions;
  size_t nnodes = problem.n;
  size_t dim = 2;
  Eigen::MatrixXf pos = Eigen::Map<const Eigen::MatrixXf>(_pos.data(), dim, nnodes);

  const double k = std::sqrt(1.0 / nnodes);
  double t = (pos.rowwise().maxCoeff() - pos.rowwise().minCoeff()).maxCoeff() * 0.1;
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
    if (delta_pos.norm() / nnodes < threshold) break;
    if (iteration % 10 == 0)
      positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), pos.size()));
  }

  positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), pos.size()));

  return positions;
}
