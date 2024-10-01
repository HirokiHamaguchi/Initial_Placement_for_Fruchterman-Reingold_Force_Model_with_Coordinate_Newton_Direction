#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "util/hex.hpp"
#include "util/problem.hpp"

void addHist(int loopCnt, const int LOOP_CNT, const bool measureTime, Grid& grid,
             std::vector<Grid>& grids) {
  if (measureTime) return;
  // if (loopCnt % problem.n == 0) grids.push_back(grid);
  if (loopCnt % (LOOP_CNT / 20) == 0) grids.push_back(grid);
}

std::pair<double, double> computeDxDy(double gx, double gy, double hxx, double hxy,
                                      double hyy) {
  double det = hxx * hyy - hxy * hxy;
  assert(det >= -1e-9);
  double inv_det = 1.0 / det;
  double newton_x = inv_det * (-hyy * gx + hxy * gy);
  double newton_y = inv_det * (-hxx * gy + hxy * gx);
  return {newton_x, newton_y};
}

std::vector<Eigen::VectorXf> solve_init(const Problem& problem,
                                        const bool measureTime) {
  auto t0 = std::chrono::high_resolution_clock::now();

  Grid grid(problem.n, problem.k);
  std::vector<Grid> grids = {grid};

  std::mt19937 gen(0);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int LOOP_CNT = problem.n * 100;
  for (int loopCnt = 0; loopCnt < LOOP_CNT; loopCnt++) {
    // randomly select a vertex
    int i = dist(gen);

    // calculate gradient and Hessian
    double gx = 0.0, gy = 0.0;
    double hxx = 0.0, hxy = 0.0, hyy = 0.0;
    for (auto [j, w] : problem.adj[i]) {
      int dq = grid.points[i].q - grid.points[j].q;
      int dr = grid.points[i].r - grid.points[j].r;
      grid.calc_grad_hess(dq, dr, w, gx, gy, hxx, hxy, hyy);
    }

    // compute Newton's direction (Hess^{-1} @ (-grad))
    auto [dx, dy] = computeDxDy(gx, gy, hxx, hxy, hyy);

    // move the vertex
    auto [x, y] = grid.hex2xy(i);
    Hex new_v = grid.xy2hex(x + dx, y + dy);
    if (grid.points[i] == new_v) {
      addHist(loopCnt, LOOP_CNT, measureTime, grid, grids);
      continue;
    }

    std::vector<Hex> path;
    for (auto& hex : grid.linedraw(grid.points[i], new_v))
      if (grid.isInside(hex)) path.push_back(hex);
    assert(std::all_of(path.begin(), path.end(),
                       [&](const Hex& hex) { return grid.isInside(hex); }));

    for (int j = 0; j < int(path.size()) - 1; ++j) {
      int& curr = grid.array[path[j].q][path[j].r];
      int& next = grid.array[path[j + 1].q][path[j + 1].r];
      if (next != -1)
        std::swap(grid.points[curr], grid.points[next]);
      else
        grid.points[i] = path[j + 1];
      std::swap(curr, next);
    }

    addHist(loopCnt, LOOP_CNT, measureTime, grid, grids);
  }
  grids.push_back(grid);

  // double kForHex = grid.computeKForHex(problem);
  std::vector<Eigen::VectorXf> positions;
  for (auto& grid : grids) {
    // grid.setKForHex(kForHex);
    positions.push_back(grid.toPosition());
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
                   1000.0
            << " seconds\n";

  return positions;
}
