#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "../util/computeDxDy.hpp"
#include "../util/grid.hpp"
#include "../util/hex.hpp"
#include "../util/problem.hpp"
#include "../util/timer.hpp"

void addVis(const Grid& grid, std::vector<Eigen::VectorXf>& positions, int it,
            bool measureTime) {
  if (measureTime) return;
  if (it % grid.n != 0) return;
  positions.push_back(grid.toPosition());
}

std::vector<Eigen::VectorXf> solve_init(const Problem& problem, const bool measureTime,
                                        const int seed,
                                        std::vector<Eigen::VectorXf>& positions,
                                        std::vector<std::pair<double, double>>& hist,
                                        Timer& timer) {
  assert(positions.empty());
  timer.start();

  Grid grid(problem.n, problem.k, seed);

  // initialize random number generator
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> distVertex(0, problem.n - 1);
  std::uniform_real_distribution<double> distHexR(0, 2.0 * grid.k);
  std::uniform_real_distribution<double> distHexTheta(0, 2 * M_PI);
  std::uniform_real_distribution<double> distSA(0.0, 1.0);

  // simulated annealing parameters
  const int ITERATIONS = 0.5 * problem.n * (problem.n * problem.n / problem.m);
  const double T0 = +2;
  const double T1 = -4;

  for (int it = 0; it < ITERATIONS; it++) {
    // randomly select a vertex
    const size_t i = distVertex(gen);
    const Hex hexI = grid.points[i];

    // calculate gradient and Hessian
    double gx = 0.0, gy = 0.0;
    double hxx = 0.0, hxy = 0.0, hyy = 0.0;
    for (auto [j, w] : problem.adj[i]) {
      int dq = hexI.q - grid.points[j].q;
      int dr = hexI.r - grid.points[j].r;
      grid.calc_grad_hess(dq, dr, w, gx, gy, hxx, hxy, hyy);
    }

    // if local minimum, continue
    if (gx * gx + gy * gy < 1e-9) {
      addVis(grid, positions, it, measureTime);
      continue;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad))
    const auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy, problem.k);

    // select random neighbor of (x + dx, y + dy)
    const auto [x, y] = grid.hex2xy(hexI.q, hexI.r);
    double r = distHexR(gen), theta = distHexTheta(gen);
    double dxr = r * std::cos(theta), dyr = r * std::sin(theta);
    const Hex hexJ = grid.xy2hex(x + dx + dxr, y + dy + dyr);
    if (!grid.isInside(hexJ)) {
      addVis(grid, positions, it, measureTime);
      continue;
    }

    // compute (postScore - preScore)
    double diffScore = grid.calcScoreDiff(problem, hexI, hexJ);

    // if accept by Simulated Annealing, continue
    double T = std::pow(10, T0 + (T1 - T0) * it / ITERATIONS);
    if (diffScore <= 0 || distSA(gen) < std::exp(-diffScore / T)) {
      // swap position
      grid.swap(i, hexI, hexJ);
      addVis(grid, positions, it, measureTime);
      continue;
    }

    // add copy of current position for visualization
    addVis(grid, positions, it, measureTime);

    // assert(grid.isCorrectState());  // ! For debug
  }

  auto finalPos = grid.toPosition();
  problem.optimalScaling(finalPos);
  positions.push_back(finalPos);
  assert(!measureTime || positions.size() == 1);

  timer.stop();
  hist.emplace_back(problem.calcScore(finalPos, true), timer.sec());

  return positions;
}
