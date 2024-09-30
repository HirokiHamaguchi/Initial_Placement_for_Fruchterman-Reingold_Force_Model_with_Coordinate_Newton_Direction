#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>  // for sqrt function
#include <iostream>
#include <random>

#include "util/hex.hpp"
#include "util/problem.hpp"

std::pair<double, double> computeEigenvalues(double h_xx, double h_xy, double h_yy) {
  // The characteristic equation is of the form: λ^2 - trace * λ + determinant = 0
  // Solving this quadratic equation: λ = (trace ± sqrt(trace^2 - 4 * determinant)) / 2

  double trace = h_xx + h_yy;
  double determinant = h_xx * h_yy - h_xy * h_xy;
  double discriminant = trace * trace - 4 * determinant;
  assert(discriminant >= -1e-9);  // The symmetric matrix should have real eigenvalues

  double lambda1 = (trace + std::sqrt(discriminant + 1e-9)) / 2.0;
  double lambda2 = (trace - std::sqrt(discriminant + 1e-9)) / 2.0;

  return {lambda1, lambda2};
}

void addHist(int loopCnt, const int LOOP_CNT, const bool measureTime,
             const Problem& problem, const Grid& grid,
             std::vector<Eigen::VectorXf>& positions) {
  if (measureTime) return;
  if (loopCnt % problem.n == 0) {
    dbg(grid.calcScore(problem, false));
    positions.push_back(grid.toPosition());
  }
}

std::pair<double, double> computeDxDy(double grad_x, double grad_y, double hess_xx,
                                      double hess_xy, double hess_yy) {
  double det = hess_xx * hess_yy - hess_xy * hess_xy;
  assert(det >= -1e-9);

  auto [lambda1, lambda2] = computeEigenvalues(hess_xx, hess_xy, hess_yy);
  assert(lambda1 >= -1e-9 && lambda2 >= -1e-9);

  double inv_det = 1.0 / det;
  double newton_x = inv_det * (-hess_yy * grad_x + hess_xy * grad_y);
  double newton_y = inv_det * (-hess_xx * grad_y + hess_xy * grad_x);
  return {newton_x, newton_y};
}

std::vector<Eigen::VectorXf> solve_init(const Problem& problem,
                                        const bool measureTime) {
  auto t0 = std::chrono::high_resolution_clock::now();

  Grid grid(problem.n, problem.k);
  std::vector<Eigen::VectorXf> positions;

  std::mt19937 gen(0);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int LOOP_CNT = problem.n * 300;
  for (int loopCnt = 0; loopCnt < LOOP_CNT; loopCnt++) {
    // randomly select a vertex
    int i = dist(gen);

    // calculate gradient and Hessian
    double grad_x = 0.0, grad_y = 0.0;
    double hess_xx = 0.0, hess_xy = 0.0, hess_yy = 0.0;
    for (auto [j, w] : problem.adj[i]) {
      int dq = grid.points[i].q - grid.points[j].q;
      int dr = grid.points[i].r - grid.points[j].r;
      auto [gx, gy, hxx, hxy, hyy] = grid.calc_grad_hess(dq, dr, w);
      grad_x += gx;
      grad_y += gy;
      hess_xx += hxx;
      hess_xy += hxy;
      hess_yy += hyy;
    }

    // compute Newton's direction (Hess^{-1} @ (-grad))
    double dx, dy;
    std::tie(dx, dy) = computeDxDy(grad_x, grad_y, hess_xx, hess_xy, hess_yy);

    // move the vertex
    auto [x, y] = grid.points[i].hex2xy(problem.k);
    Hex new_v = Hex::xy2hex(x + dx, y + dy, problem.k);
    if (grid.points[i] == new_v) {
      addHist(loopCnt, LOOP_CNT, measureTime, problem, grid, positions);
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

    addHist(loopCnt, LOOP_CNT, measureTime, problem, grid, positions);
  }
  positions.push_back(grid.toPosition());

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
                   1000.0
            << " seconds\n";

  return positions;
}
