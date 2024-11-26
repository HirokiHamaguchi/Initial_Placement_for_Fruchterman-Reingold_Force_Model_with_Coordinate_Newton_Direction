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

void solve_init(const Problem& problem, const bool measureTime, const int seed,
                std::vector<Eigen::VectorXf>& positions, Timer& timer) {
  timer.start();

  Grid grid(problem.n, problem.k, seed);

  // initialize random number generator
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> distVertex(0, problem.n - 1);
  std::uniform_real_distribution<double> distHexR(0, 1);
  std::uniform_real_distribution<double> distHexTheta(0, 2 * M_PI);

  // simulated annealing parameters
  const int ITERATIONS = 2 * problem.n * (problem.n * problem.n / problem.m);
  const double T0 = +1.5;
  const double T1 = +0;

  for (int it = 0; it < ITERATIONS; it++) {
    // randomly select a vertex
    const size_t i = distVertex(gen);
    const Hex hexI = grid.points[i];

    // calculate gradient and Hessian
    float gx = 0.0, gy = 0.0;
    float hxx = 0.0, hxy = 0.0, hyy = 0.0;
    for (auto [j, w] : problem.adj[i]) {
      int dq = hexI.q - grid.points[j].q;
      int dr = hexI.r - grid.points[j].r;
      grid.calc_grad_hess(dq, dr, w, gx, gy, hxx, hxy, hyy);
    }

    // if local minimum, continue
    if ((gx * gx + gy * gy) * grid.k * grid.k < 1e-9) {
      addVis(grid, positions, it, measureTime);
      continue;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad))
    const auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy, problem.k);

    // select random neighbor of (x + dx, y + dy) with randomness
    const auto [x, y] = grid.hex2xy(hexI.q, hexI.r);
    float T = T0 + (T1 - T0) * it / ITERATIONS;
    float r = T * distHexR(gen), theta = distHexTheta(gen);
    float dxr = r * std::cos(theta), dyr = r * std::sin(theta);
    const Hex hexJ = grid.xy2hex(x + dx + dxr, y + dy + dyr);

    // swap position
    if (grid.isInside(hexJ) && hexI != hexJ) grid.swap(i, hexI, hexJ);

    addVis(grid, positions, it, measureTime);

    // assert(grid.isCorrectState());  // ! For debug
  }

  auto finalPos = grid.toPosition();
  problem.optimalScaling(finalPos);
  positions.push_back(finalPos);

  timer.stop();
}
