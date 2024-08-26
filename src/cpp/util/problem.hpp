#pragma once

#include <iostream>
#include <utility>
#include <vector>

struct Problem {
  size_t n;                  // number of vertices
  size_t m;                  // number of edges
  double k;                  // constant for score calculation
  std::vector<size_t> row;   // edge u
  std::vector<size_t> col;   // edge v
  std::vector<double> data;  // edge weight
  size_t t;                  // number of turns

  Problem(size_t n, size_t m, double k, std::vector<size_t> row,
          std::vector<size_t> col, std::vector<double> data, size_t t)
      : n(n), m(m), k(k), row(row), col(col), data(data), t(t) {}

  // positions[turn][vertex]
  void printOutput(std::vector<std::vector<std::pair<double, double>>> positions) {
    std::cout << n << " " << m << " " << k << "\n";
    for (size_t i = 0; i < m; ++i) {
      std::cout << row[i] << " " << col[i] << " " << data[i] << "\n";
    }
    std::cout << t << "\n";
    for (const auto& turn : positions)
      for (const auto& pos : turn) std::cout << pos.first << " " << pos.second << "\n";
  }
};
