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

bool isFinished(const Problem& problem, const Grid& grid,
                std::vector<Eigen::VectorXf>& positions, size_t& fail,
                std::vector<std::pair<double, Eigen::VectorXf>>& localHist, int it,
                bool measureTime, bool isFailed) {
  // if the vertex is not moved, increment fail
  if (isFailed) {
    fail++;
    if (fail >= grid.n) {
      std::cerr << "Early stop (fail = " << fail << ")" << std::endl;
      return true;
    }
  }
  fail = 0;

  // break the loop if the score is not improved
  if (it % (50 * grid.n) == 0) {
    auto pos = grid.toPosition();

    // add for visualization
    if (!measureTime) positions.push_back(grid.toPosition());

    // add for score history
    double nowScore = problem.calcScore(pos, false);
    localHist.emplace_back(nowScore, pos);
    size_t sz = localHist.size();
    if (int(sz) >= 3 && localHist[sz - 3].first <= localHist[sz - 2].first &&
        localHist[sz - 3].first <= localHist[sz - 1].first) {
      std::cerr << "Early stop (score = " << nowScore << ")" << std::endl;
      if (!measureTime) positions.push_back(localHist[sz - 3].second);
      return true;
    }
  }

  return false;
}

std::vector<Eigen::VectorXf> solve_init(const Problem& problem, const bool measureTime,
                                        const int seed,
                                        std::vector<Eigen::VectorXf>& positions,
                                        std::vector<std::pair<double, double>>& hist,
                                        Timer& timer) {
  timer.start();

  Grid grid(problem.n, problem.k);
  assert(positions.empty());
  if (!measureTime) positions.emplace_back(grid.toPosition());

  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int ITERATIONS = 3 * problem.n * (problem.n * problem.n / problem.m);
  size_t fail = 0;
  std::vector<std::pair<double, Eigen::VectorXf>> localHist;
  for (int it = 0; it < ITERATIONS; it++) {
    // randomly select a vertex
    size_t i = dist(gen);

    // calculate gradient and Hessian
    double gx = 0.0, gy = 0.0;
    double hxx = 0.0, hxy = 0.0, hyy = 0.0;
    for (auto [j, w] : problem.adj[i]) {
      int dq = grid.points[i].q - grid.points[j].q;
      int dr = grid.points[i].r - grid.points[j].r;
      grid.calc_grad_hess(dq, dr, w, gx, gy, hxx, hxy, hyy);
    }
    if (std::abs(gx) + std::abs(gy) < 1e-9) {
      if (isFinished(problem, grid, positions, fail, localHist, it, measureTime, true))
        break;
      continue;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad)), if possible
    auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy, problem.k);

    // move the vertex
    auto [x, y] = grid.hex2xy(i);
    Hex new_v = grid.xy2hex(x + dx, y + dy);
    if (grid.points[i] == new_v) {
      if (isFinished(problem, grid, positions, fail, localHist, it, measureTime, true))
        break;
      continue;
    }

    // update along path
    grid.updateAlongPath(i, new_v);

    if (isFinished(problem, grid, positions, fail, localHist, it, measureTime, false))
      break;
  }
  auto finalPos = grid.toPosition();
  problem.optimalScaling(finalPos);
  positions.push_back(finalPos);
  assert(!measureTime || positions.size() == 1);

  timer.stop();
  hist.emplace_back(problem.calcScore(finalPos, true), timer.sec());

  return positions;
}
