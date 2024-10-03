#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "solve_init.cpp"
#include "util/hex.hpp"
#include "util/problem.hpp"

std::tuple<double, double, double, double, double> computeGradHess(
    const Problem& problem, const Eigen::MatrixXf& pos, size_t i) {
  double gx = 0.0, gy = 0.0;
  double hxx = 0.0, hxy = 0.0, hyy = 0.0;
  for (auto [j, w] : problem.adj[i]) {
    auto delta = pos.col(i) - pos.col(j);
    double dist = delta.norm();
    assert(dist > 1e-9);

    double coeff1 = w * dist / problem.k;
    double coeff2 = w / (dist * problem.k);

    gx += coeff1 * delta[0];
    gy += coeff1 * delta[1];
    hxx += coeff1 + coeff2 * delta[0] * delta[0];
    hxy += coeff2 * delta[0] * delta[1];
    hyy += coeff1 + coeff2 * delta[1] * delta[1];
  }
  for (size_t j = 0; j < problem.n; j++) {
    if (j == i) continue;
    auto delta = pos.col(i) - pos.col(j);
    double dist = delta.norm();
    assert(dist > 1e-9);

    double d2 = std::pow(dist, 2);
    double d4 = std::pow(d2, 2);
    double k2 = std::pow(problem.k, 2);
    double coeff1 = -k2 / d2;
    double coeff2 = 2 * k2 / d4;

    gx += coeff1 * delta[0];
    gy += coeff1 * delta[1];
    hxx += coeff1 + coeff2 * delta[0] * delta[0];
    hxy += coeff2 * delta[0] * delta[1];
    hyy += coeff1 + coeff2 * delta[1] * delta[1];
  }
  return {gx, gy, hxx, hxy, hyy};
}

std::vector<Eigen::VectorXf> solve_RS(const Problem& problem,
                                      const Eigen::VectorXf& _pos,
                                      const bool measureTime) {
  auto t0 = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXf pos = Eigen::Map<const Eigen::MatrixXf>(_pos.data(), 2, problem.n);
  std::mt19937 gen(0);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int LOOP_CNT = 10 * problem.n;

  std::vector<Eigen::VectorXf> positions;

  for (int loopCnt = 0; loopCnt < LOOP_CNT; loopCnt++) {
    // randomly select a vertex
    size_t i = dist(gen);

    // calculate gradient and Hessian
    auto [gx, gy, hxx, hxy, hyy] = computeGradHess(problem, pos, i);

    // compute Newton's direction (Hess^{-1} @ (-grad), if Hess is positive definite)
    auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy, problem.k);

    if (std::hypot(dx, dy) > problem.k) {
      dbg(i, gx, gy, hxx, hxy, hyy, dx, dy);
    }

    // move the vertex
    pos.col(i) += Eigen::Vector2f(dx, dy);

    if (!measureTime && loopCnt % problem.n == 0)
      positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), 2 * problem.n));
  }
  positions.push_back(Eigen::Map<Eigen::VectorXf>(pos.data(), 2 * problem.n));

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "RS Elapsed time: " << std::chrono::duration<double>(t1 - t0).count()
            << " seconds\n";

  return positions;
}
