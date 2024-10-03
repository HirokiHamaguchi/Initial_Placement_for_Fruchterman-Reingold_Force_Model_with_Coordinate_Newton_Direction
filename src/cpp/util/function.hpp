#pragma once

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "problem.hpp"

/**
 * @brief struct for L-BFGS solver (Fruchterman-Reingold layout)
 * @note operator() returns the score and gradient
 */
struct FunctionFR {
  const Problem& problem;
  FunctionFR(const Problem& problem) : problem(problem) {}

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

/**
 * @brief struct for L-BFGS solver (Kamada-Kawai layout)
 * @note operator() returns the score and gradient
 */
struct FunctionKK {
  const Problem& problem;
  Eigen::MatrixXf invdist;

  FunctionKK(const Problem& problem) : problem(problem) {
    // calculate distance between vertices
    invdist.resize(problem.n, problem.n);
    for (size_t u = 0; u < problem.n; ++u)
      invdist.col(u) = dijkstra(u).array().inverse();
    for (size_t u = 0; u < problem.n; ++u) invdist(u, u) = 1e+3;
  }

  Eigen::VectorXf dijkstra(size_t start) {
    Eigen::VectorXf dist(problem.n);
    dist.setConstant(std::numeric_limits<float>::infinity());
    dist[start] = 0.0;

    std::priority_queue<std::pair<float, size_t>, std::vector<std::pair<float, size_t>>,
                        std::greater<std::pair<float, size_t>>>
        pq;
    pq.push({0.0, start});

    while (!pq.empty()) {
      auto [d, u] = pq.top();
      pq.pop();
      if (dist[u] < d) continue;
      for (auto [v, w] : problem.adj[u]) {
        float d2 = d + w;
        if (d2 < dist[v]) {
          dist[v] = d2;
          pq.push({d2, v});
        }
      }
    }

    return dist;
  }

  float operator()(const Eigen::VectorXf& x, Eigen::VectorXf& grad) {
    grad.setZero();

    int dim = 2;
    int nNodes = problem.n;
    double meanWeight = 1e-3;

    Eigen::MatrixXf pos =
        Eigen::Map<const Eigen::MatrixXf>(x.data(), dim, nNodes).transpose();

    double cost = 0.0;
    for (int i = 0; i < nNodes; ++i) {
      Eigen::MatrixXf delta = pos.row(i).replicate(nNodes, 1) - pos;  // n times 2
      Eigen::VectorXf nodeSep = delta.rowwise().norm();               // n times 1
      nodeSep[i] += 1e-3;
      Eigen::MatrixXf direction =
          delta.array() / nodeSep.replicate(1, dim).array();  // n times 2
      nodeSep[i] -= 1e-3;
      Eigen::VectorXf offset = nodeSep.cwiseProduct(invdist.col(i)) -
                               Eigen::VectorXf::Ones(nNodes);  // n times 1
      offset[i] = 0;
      cost += 0.5 * offset.array().square().sum();
      auto offsetInvdist = offset.cwiseProduct(invdist.col(i));  // n times 1
      Eigen::MatrixXf offsetInvdistDir = offsetInvdist.replicate(1, dim)
                                             .cwiseProduct(direction)
                                             .transpose();  // 2 times n
      assert(offsetInvdistDir.rows() == 2);
      assert(offsetInvdistDir.cols() == nNodes);
      grad.segment<2>(2 * i) += direction.transpose() * offsetInvdist;
      grad -= Eigen::Map<const Eigen::VectorXf>(offsetInvdistDir.data(), 2 * nNodes);
    }

    // Additional parabolic term to encourage mean position to be near origin:
    Eigen::VectorXf sumPos = pos.colwise().sum();
    cost += 0.5 * meanWeight * sumPos.squaredNorm();
    Eigen::MatrixXf gradMat = meanWeight * sumPos.replicate(1, nNodes);
    grad += Eigen::Map<const Eigen::VectorXf>(gradMat.data(), 2 * nNodes);

    return cost;
  }
};
