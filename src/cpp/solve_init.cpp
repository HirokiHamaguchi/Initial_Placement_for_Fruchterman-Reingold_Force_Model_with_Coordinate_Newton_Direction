#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <random>

#include "util/hex.hpp"
#include "util/problem.hpp"

std::vector<Eigen::VectorXf> solve_init(const Problem& problem) {
  auto t0 = std::chrono::high_resolution_clock::now();

  Grid grid(problem.n, problem.k);
  std::vector<Eigen::VectorXf> positions;

  std::mt19937 gen(0);
  std::uniform_int_distribution<int> dist(0, problem.n - 1);
  const int LOOP_CNT = problem.n * 1000;
  for (int loopCnt = 0; loopCnt < LOOP_CNT; loopCnt++) {
    // randomly select a vertex
    int i = dist(gen);

    // compute Newton's direction of the vertex
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

    // calculate Newton's direction (Hess^{-1} @ (-grad))
    double det = hess_xx * hess_yy - hess_xy * hess_xy;
    assert(det > 1e-9);
    double inv_det = 1.0 / det;
    double newton_x = inv_det * (-hess_yy * grad_x + hess_xy * grad_y);
    double newton_y = inv_det * (-hess_xx * grad_y + hess_xy * grad_x);
    if (newton_x == 0.0 && newton_y == 0.0) continue;

    // move the vertex
    auto [x, y] = grid.points[i].hex2xy(problem.k);
    Hex new_v = grid.points[i];
    while (grid.points[i] == new_v) {
      new_v = Hex::xy2hex(x + newton_x, y + newton_y, problem.k);
      newton_x *= 1.5;
      newton_y *= 1.5;
    }
    assert(grid.isInside(new_v));
    std::vector<Hex> path = grid.linedraw(grid.points[i], new_v);

    // assert(path.size() >= 2);
    // if (path.size() == 2) {
    //   std::vector<int> vertices;
    //   for (auto& hex : path) {
    //     int v = grid.array[hex.q][hex.r];
    //     if (v != -1) vertices.push_back(v);
    //   }
    //   assert(std::count(vertices.begin(), vertices.end(), i) == 1);

    // double scoreBefore = 0, scoreAfter = 0;
    // for (int v : vertices) {
    //   for (auto [u, w] : problem.adj[v]) {
    //     auto [dx, dy] = (grid.points[v] - grid.points[u]).hex2xy(problem.k);
    //     double d = std::hypot(dx, dy);
    //     scoreBefore += w * std::pow(d, 3) / (3.0 * problem.k);
    //   }
    // }

    // int& curr = grid.array[path[0].q][path[0].r];
    // int& next = grid.array[path[1].q][path[1].r];
    // if (next != -1)
    //   std::swap(grid.points[curr], grid.points[next]);
    // else
    //   grid.points[i] = path[1];
    // std::swap(curr, next);

    // for (int v : vertices) {
    //   for (auto [u, w] : problem.adj[v]) {
    //     auto [dx, dy] = (grid.points[v] - grid.points[u]).hex2xy(problem.k);
    //     double d = std::hypot(dx, dy);
    //     scoreAfter += w * std::pow(d, 3) / (3.0 * problem.k);
    //   }
    // }

    // if (scoreAfter > scoreBefore) {
    //   // dbg("rollback");
    //   std::swap(curr, next);
    //   if (next != -1)
    //     std::swap(grid.points[curr], grid.points[next]);
    //   else
    //     grid.points[i] = path[0];
    // } else {
    //   // dbg("move");
    // }
    // } else {
    // move vertex along path
    for (int j = 0; j < int(path.size()) - 1; ++j) {
      int& curr = grid.array[path[j].q][path[j].r];
      int& next = grid.array[path[j + 1].q][path[j + 1].r];
      if (next != -1)
        std::swap(grid.points[curr], grid.points[next]);
      else
        grid.points[i] = path[j + 1];
      std::swap(curr, next);
    }
    // }

    // add hist
    if (loopCnt % (LOOP_CNT / 30) == 0) {
      dbg(grid.calcScore(problem, false));
      assert(grid.isCorrectState());
      Eigen::VectorXf position(2 * problem.n);
      for (int i = 0; i < int(problem.n); ++i)
        std::tie(position[2 * i], position[2 * i + 1]) =
            grid.points[i].hex2xy(problem.k);
      positions.push_back(position);
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
                   1000.0
            << " seconds\n";

  return positions;
}
