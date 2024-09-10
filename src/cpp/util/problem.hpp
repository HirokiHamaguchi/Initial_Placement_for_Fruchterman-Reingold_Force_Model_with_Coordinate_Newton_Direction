#pragma once

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "dbg.h"

struct Problem {
  size_t n;                         // number of vertices
  size_t m;                         // number of edges
  double k;                         // constant for score calculation
  std::vector<size_t> row;          // edge u
  std::vector<size_t> col;          // edge v
  std::vector<double> data;         // edge weight
  std::string matrixName = "test";  // matrix name

  Problem(const std::string matrixName) : matrixName(matrixName) {
    std::string curPath = std::filesystem::current_path().string();
    assert(curPath.substr(curPath.size() - 8, 8) == "\\src\\cpp");
    std::string path =
        curPath.substr(0, curPath.size() - 8) + "\\data\\" + matrixName + ".mtx";

    std::ifstream file(path);
    if (!file.is_open()) {
      std::cerr << "Error: file not found\n";
      exit(1);
    }

    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '%') continue;
      break;
    }

    size_t _n, _m;
    std::istringstream iss(line);
    iss >> n >> _n >> _m;
    assert(n == _n);

    std::map<std::pair<size_t, size_t>, double> edges;
    for (size_t i = 0; i < _m; ++i) {
      std::getline(file, line);
      std::istringstream iss(line);
      int r, c;
      double w;
      iss >> r >> c >> w;
      r--;
      c--;
      if (r == c) continue;
      edges[std::minmax(r, c)] += std::abs(w);
    }

    m = edges.size();
    row.reserve(m);
    col.reserve(m);
    data.reserve(m);

    for (const auto& [rc, w] : edges) {
      row.push_back(rc.first);
      col.push_back(rc.second);
      data.push_back(w);
    }

    file.close();
    k = 1 / std::sqrt(n);
  }

  Problem(size_t n, size_t m, double k, std::vector<size_t> row,
          std::vector<size_t> col, std::vector<double> data)
      : n(n), m(m), k(k), row(row), col(col), data(data) {}

  double calcScore(const std::vector<std::pair<double, double>>& pos,
                   bool addRepulsive = true) const {
    double score = 0.0;

    if (addRepulsive) {
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
          double x1 = pos[i].first;
          double y1 = pos[i].second;
          double x2 = pos[j].first;
          double y2 = pos[j].second;
          double d = std::hypot(x1 - x2, y1 - y2);
          if (d < 1e-9) {
            dbg(i, j, x1, y1, x2, y2, d);
            assert(false);
          }
          score -= std::pow(k, 2) * std::log(d);
        }
      }
    }

    for (size_t i = 0; i < m; ++i) {
      size_t u = row[i];
      size_t v = col[i];
      double a = data[i];
      double x1 = pos[u].first;
      double y1 = pos[u].second;
      double x2 = pos[v].first;
      double y2 = pos[v].second;
      double d = std::hypot(x1 - x2, y1 - y2);
      score += a * std::pow(d, 3) / (3.0 * k);
    }

    return score;
  }

  // positions[turn][vertex]
  void printOutput(std::vector<std::vector<std::pair<double, double>>> positions) {
    std::string curPath = std::filesystem::current_path().string();
    assert(curPath.substr(curPath.size() - 8, 8) == "\\src\\cpp");
    std::string path =
        curPath.substr(0, curPath.size() - 8) + "\\out\\" + matrixName + ".out";
    std::ofstream file(path);
    if (!file.is_open()) {
      dbg(path);
      std::cerr << "Error: file not found\n";
      exit(1);
    }

    file << n << " " << m << " " << k << "\n";
    for (size_t i = 0; i < m; ++i) {
      file << row[i] << " " << col[i] << " " << data[i] << "\n";
    }
    file << positions.size() << "\n";
    for (const auto& turn : positions)
      for (const auto& pos : turn) file << pos.first << " " << pos.second << "\n";
    file.close();

    std::cout << "Output path:" << path << "\n";
  }
};
