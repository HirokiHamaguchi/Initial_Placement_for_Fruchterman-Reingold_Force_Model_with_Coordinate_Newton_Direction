#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include "dbg.h"
#include "openFile.hpp"
#include "forceModel.hpp"
#include "frModel.hpp"
#include "hcModel.hpp"
#include "eadesModel.hpp"

enum class ForceModelType
{
  FR,
  HC,
  Eades
};
struct Problem
{
  size_t n;                                                // number of vertices
  size_t m;                                                // number of edges
  double k;                                                // constant for score calculation
  std::vector<size_t> row;                                 // edge u
  std::vector<size_t> col;                                 // edge v
  std::vector<double> data;                                // edge weight
  std::vector<std::vector<std::pair<size_t, double>>> adj; // adjacency list
  std::string matrixName = "test";                         // matrix name
  std::unique_ptr<ForceModel> model;                       // force model implementation
  ForceModelType modelType;                                // which model in use

  Problem() : n(0), m(0), k(0.0), modelType(ForceModelType::FR) {}
  Problem(size_t n, size_t m, double k, std::vector<size_t> &row,
          std::vector<size_t> &col, std::vector<double> &data,
          ForceModelType _modelType = ForceModelType::FR)
      : n(n), m(m), k(k), row(row), col(col), data(data), modelType(_modelType)
  {
    assert(row.size() == m);
    assert(col.size() == m);
    assert(data.size() == m);
    makeAdj();
    createModel();
    assert(isConnected());
  }

  Problem(const std::string matrixName,
          ForceModelType _modelType = ForceModelType::FR,
          bool is1 = false)
      : matrixName(matrixName), modelType(_modelType)
  {
    // if is1, then all edge weights are set to 1

    std::string curPath = std::filesystem::current_path().string();
    size_t folderPathSize = curPath.find("Initial_Placement_for_FR");
    assert(folderPathSize != std::string::npos);
    curPath = curPath.substr(
        0, folderPathSize + std::string("Initial_Placement_for_FR").size());
    std::string path = curPath + "/data/" + matrixName + ".mtx";
    assert(std::filesystem::exists(path));

    std::ifstream file(path);
    if (!file.is_open())
    {
      std::cerr << "Error: file not found\n";
      exit(1);
    }

    std::string line;
    std::getline(file, line);
    assert(std::count(line.begin(), line.end(), ' ') == 4);
    std::vector<std::string> tokens;
    std::istringstream iss0(line);
    for (std::string token; iss0 >> token;)
      tokens.push_back(token);
    assert(tokens[0] == "%MatrixMarket" || tokens[0] == "%%MatrixMarket");
    assert(tokens[1] == "matrix");
    assert(tokens[2] == "coordinate");
    std::string field = tokens[3];
    assert(field == "pattern" || field == "real" || field == "integer");

    while (std::getline(file, line))
    {
      if (line[0] == '%')
        continue;
      break;
    }

    size_t _n, _m;
    std::istringstream iss(line);
    iss >> n >> _n >> _m;
    assert(n == _n);

    std::map<std::pair<size_t, size_t>, double> edges;
    for (size_t i = 0; i < _m; ++i)
    {
      std::getline(file, line);
      std::istringstream iss(line);
      if (field == "pattern")
      {
        int r, c;
        iss >> r >> c;
        r--;
        c--;
        if (r == c)
          continue;
        edges[std::minmax(r, c)] += 1;
      }
      else
      {
        int r, c;
        double w;
        iss >> r >> c >> w;
        r--;
        c--;
        if (r == c)
          continue;
        edges[std::minmax(r, c)] += std::abs(w);
      }
    }
    // set all edge weights to 1
    if (is1)
      for (auto &[_, w] : edges)
        w = 1;

    m = edges.size();
    row.reserve(m);
    col.reserve(m);
    data.reserve(m);

    for (const auto &[rc, w] : edges)
    {
      row.push_back(rc.first);
      col.push_back(rc.second);
      data.push_back(w);
    }

    file.close();

    k = 1.0 / std::sqrt(n);
    makeAdj();
    createModel();
    assert(isConnected());
  }

  void createModel()
  {
    switch (modelType)
    {
    case ForceModelType::FR:
      model = std::make_unique<FRModel>(n, m, k, row, col, data);
      break;
    case ForceModelType::HC:
      model = std::make_unique<HCModel>(n, m, row, col, data);
      break;
    case ForceModelType::Eades:
      model = std::make_unique<EadesModel>(n, m, row, col, data);
      break;
    }
  }

  void makeAdj()
  {
    adj.resize(n);
    for (size_t i = 0; i < m; ++i)
    {
      adj[row[i]].emplace_back(col[i], data[i]);
      adj[col[i]].emplace_back(row[i], data[i]);
    }
  }

  bool isConnected() const
  {
    std::vector<bool> visited(n, false);
    auto dfs = [&](const auto &f, size_t u) -> void
    {
      visited[u] = true;
      for (auto [v, _] : adj[u])
        if (!visited[v])
          f(f, v);
    };
    dfs(dfs, 0);

    return std::all_of(visited.begin(), visited.end(), [](bool v)
                       { return v; });
  }

  void printOutput(const std::vector<Eigen::VectorXd> &positions,
                   std::string _path) const
  {
    auto [path, file] = openFile(_path);

    file << n << " " << m << " " << k << "\n";
    for (size_t i = 0; i < m; ++i)
    {
      file << row[i] << " " << col[i] << " " << data[i] << "\n";
    }
    file << positions.size() << "\n";
    for (auto &position : positions)
      for (size_t i = 0; i < n; ++i)
        file << position[2 * i] << " " << position[2 * i + 1] << "\n";
    file.close();

    std::cout << "Output path: " << path << "\n";
  }

  double calcScore(const Eigen::VectorXd &position,
                   bool includeRepulsive = true) const
  {
    return model->calcScore(position, includeRepulsive);
  }

  double calc_score_and_grad(const Eigen::VectorXd &x, Eigen::VectorXd &grad) const
  {
    return model->calc_score_and_grad(x, grad);
  }

  void calc_grad_hess(double dist, double dx, double dy, double w, double &gx, double &gy, double &hxx,
                      double &hxy, double &hyy) const
  {
    model->calc_grad_hess(dist, dx, dy, w, gx, gy, hxx, hxy, hyy);
  }

  double optimalScaling(Eigen::VectorXd &position) const
  {
    return model->optimalScaling(position);
  }
};
