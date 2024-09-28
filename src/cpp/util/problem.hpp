#pragma once

#include <Eigen/Core>
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
  size_t n;                  // number of vertices
  size_t m;                  // number of edges
  double k;                  // constant for score calculation
  std::vector<size_t> row;   // edge u
  std::vector<size_t> col;   // edge v
  std::vector<double> data;  // edge weight
  std::vector<std::vector<std::pair<size_t, double>>> adj;  // adjacency list
  std::string matrixName = "test";                          // matrix name

  Problem(size_t n, size_t m, double k, std::vector<size_t> row,
          std::vector<size_t> col, std::vector<double> data)
      : n(n), m(m), k(k), row(row), col(col), data(data) {
    assert(row.size() == m);
    assert(col.size() == m);
    assert(data.size() == m);
    makeAdj();
  }

  Problem(const std::string matrixName) : matrixName(matrixName) {
    std::string curPath = std::filesystem::current_path().string();
    assert(curPath.substr(curPath.size() - 8, 8) == "/src/cpp");
    std::string path =
        curPath.substr(0, curPath.size() - 8) + "/data/" + matrixName + ".mtx";

    std::ifstream file(path);
    if (!file.is_open()) {
      std::cerr << "Error: file not found\n";
      exit(1);
    }

    std::string line;
    std::getline(file, line);
    assert(std::count(line.begin(), line.end(), ' ') == 4);
    std::vector<std::string> tokens;
    std::istringstream iss0(line);
    for (std::string token; iss0 >> token;) tokens.push_back(token);
    assert(tokens[0] == "%MatrixMarket" || tokens[0] == "%%MatrixMarket");
    assert(tokens[1] == "matrix");
    assert(tokens[2] == "coordinate");
    std::string field = tokens[3];
    assert(field == "pattern" || field == "real" || field == "integer");

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
      if (field == "pattern") {
        int r, c;
        iss >> r >> c;
        r--;
        c--;
        if (r == c) continue;
        edges[std::minmax(r, c)] += 1;
      } else {
        int r, c;
        double w;
        iss >> r >> c >> w;
        r--;
        c--;
        if (r == c) continue;
        edges[std::minmax(r, c)] += std::abs(w);
      }
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
    makeAdj();
    assert(isConnected());
  }

  void makeAdj() {
    adj.resize(n);
    for (size_t i = 0; i < m; ++i) {
      adj[row[i]].emplace_back(col[i], data[i]);
      adj[col[i]].emplace_back(row[i], data[i]);
    }
  }

  bool isConnected() const {
    std::vector<bool> visited(n, false);
    auto dfs = [&](const auto& f, size_t u) -> void {
      visited[u] = true;
      for (auto [v, _] : adj[u])
        if (!visited[v]) f(f, v);
    };
    dfs(dfs, 0);

    return std::all_of(visited.begin(), visited.end(), [](bool v) { return v; });
  }

  // positions[turn][vertex]
  void printOutput(const std::vector<Eigen::VectorXf>& positions) {
    std::string curPath = std::filesystem::current_path().string();
    assert(curPath.substr(curPath.size() - 8, 8) == "/src/cpp");
    std::string path =
        curPath.substr(0, curPath.size() - 8) + "/out/" + matrixName + ".out";
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
    for (auto& position : positions)
      for (size_t i = 0; i < n; ++i)
        file << position[2 * i] << " " << position[2 * i + 1] << "\n";
    file.close();

    std::cout << "Output path:" << path << "\n";
  }
};
