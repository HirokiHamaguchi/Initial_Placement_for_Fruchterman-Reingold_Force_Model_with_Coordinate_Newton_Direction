#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <random>

#include "util/hex.hpp"
#include "util/problem.hpp"

Eigen::VectorXf solve_init(const Problem& problem) {
  auto t0 = std::chrono::high_resolution_clock::now();

  int n = problem.n;
  double k = problem.k;
  std::vector<Hex> points;
  int r = 0;
  while (int(points.size()) < n) {
    for (int q = -(r / 2); q < int(1 / k) - (r / 2); ++q) {
      Hex h(q, r, -q - r);
      points.push_back(h);
      if (int(points.size()) == n) break;
    }
    r++;
  }
  std::random_shuffle(points.begin(), points.end());

  std::mt19937 gen(0);
  std::uniform_int_distribution<int> dist(0, n - 1);
  for (int _ = 0; _ < 1000; _++) {
    // randomly select a vertex
    int v = dist(gen);
    dbg(v);

    // compute Newton's direction of the vertex
    // todo
    int x, y;
    dbg(x, y);

    // move the vertex
    Hex new_v = Hex::xy2hex(x, y, k);
    dbg(new_v);
    std::vector<Hex> path = Hex::linedraw(points[v], new_v);
    dbg(path);

    // calculate the score diff (line search with backtracking)
    while (path.size()) {
      // if failed, shrink the path to the half
      path.resize(path.size() / 2);
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::seconds>(t1 - t0).count()
            << " seconds\n";

  Eigen::VectorXf position(2 * n);
  for (int i = 0; i < n; ++i) {
    std::tie(position[2 * i], position[2 * i + 1]) = points[i].xy(k);
  }
  return position;
}
