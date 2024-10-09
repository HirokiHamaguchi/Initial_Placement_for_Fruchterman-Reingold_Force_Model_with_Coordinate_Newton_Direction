#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "util/computeDxDy.hpp"
#include "util/grid.hpp"
#include "util/hex.hpp"
#include "util/problem.hpp"

bool isFinished(const Problem& problem, const Grid& grid,
                std::vector<Eigen::VectorXf>& positions, size_t& fail,
                std::vector<double>& scores, int it, bool measureTime, bool isFailed) {
  // add for visualization
  if (!measureTime && it % (20 * grid.n) == 0) positions.push_back(grid.toPosition());

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
  if (it % (20 * grid.n) == 0) {
    double nowScore = problem.calcScore(grid.toPosition(), false);
    scores.push_back(nowScore);
    // もし、scoresが3つ以上あり、かつ、最小値がscores[-3]である場合、終了
    size_t sz = scores.size();
    if (int(sz) >= 3 && scores[sz - 3] <= scores[sz - 2] &&
        scores[sz - 3] <= scores[sz - 1]) {
      std::cerr << "Early stop (score = " << nowScore << ")" << std::endl;
      // return true;
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
  positions.emplace_back(grid.toPosition());

  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int ITERATIONS = problem.n * 300;
  size_t fail = 0;
  std::vector<double> scores;
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
      if (isFinished(problem, grid, positions, fail, scores, it, measureTime, true))
        break;
      continue;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad)), if possible
    auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy, problem.k);

    // move the vertex
    auto [x, y] = grid.hex2xy(i);
    Hex new_v = grid.xy2hex(x + dx, y + dy);
    if (grid.points[i] == new_v) {
      if (isFinished(problem, grid, positions, fail, scores, it, measureTime, true))
        break;
      continue;
    }

    // update along path
    grid.updateAlongPath(i, new_v);

    if (isFinished(problem, grid, positions, fail, scores, it, measureTime, false))
      break;
  }
  positions.push_back(grid.toPosition());

  timer.stop();
  hist.emplace_back(problem.calcScore(grid.toPosition(), false), timer.sec());

  return positions;
}
