#pragma once

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "hex.hpp"
#include "problem.hpp"

struct Grid {
  size_t n;        // number of vertices
  size_t n2;       // length of the side of the hexagon
  size_t arraySz;  // size of the array
  double k;
  std::vector<Hex> points;
  std::vector<int> array;

 public:
  Grid(int n, double k) : n(n), n2(0), k(k), k2(std::pow(k, 2)) {
    size_t hexSize = 2 * n;
    while (3 * n2 * n2 + 3 * n2 + 1 < hexSize) n2++;
    arraySz = 2 * n2 + 1;
    for (int r = 0; r < int(arraySz); ++r) {
      for (int q = 0; q < int(arraySz); ++q) {
        if (r + q < int(n2) || r + q > int(3 * n2)) continue;
        points.emplace_back(q, r);
      }
    }
    assert(points.size() >= hexSize);

    std::mt19937 g(0);
    std::shuffle(points.begin(), points.end(), g);
    points.resize(n);

    array.resize(arraySz * arraySz, -1);
    for (size_t i = 0; i < points.size(); ++i)
      array[to1DIndex(points[i].q, points[i].r)] = i;
  }

  inline std::pair<float, float> hex2xy(double q, double r) const {
    return {k * (q + r / 2.0), k * (r * std::sqrt(3) / 2.0)};
  }
  inline std::pair<float, float> hex2xy(size_t i) const {
    return hex2xy(points[i].q, points[i].r);
  }
  Hex xy2hex(float x, float y) {
    float r = y * 2.0 / (k * std::sqrt(3));
    float q = x / k - r / 2.0;
    return Hex::round(q, r, -q - r);
  }

  void calc_grad_hess(int dq, int dr, double w, double& gx, double& gy, double& hxx,
                      double& hxy, double& hyy) const {
    auto delta = hex2xy(dq, dr);
    double dist = std::hypot(delta.first, delta.second);
    assert(dist > 1e-9);

    // * Only use attractive force
    double coeff1 = w * dist / k;
    double coeff2 = w / (dist * k);

    // // Method 2: use both attractive and repulsive forces (deprecated)
    // // double d2 = std::pow(dist, 2);
    // // double d4 = std::pow(d2, 2);
    // // double coeff1 = w * dist / k - k2 / d2;
    // // double coeff2 = w / (dist * k) + 2 * k2 / d4;

    gx += coeff1 * delta.first;
    gy += coeff1 * delta.second;
    hxx += coeff1 + coeff2 * delta.first * delta.first;
    hxy += coeff2 * delta.first * delta.second;
    hyy += coeff1 + coeff2 * delta.second * delta.second;
  }

  void updateAlongPath(int i, const Hex& new_v) {
    std::vector<Hex> path = linedraw(points[i], new_v);
    if (std::any_of(path.begin(), path.end(),
                    [&](const Hex& hex) { return !isInside(hex); })) {
      std::vector<Hex> newPath;
      for (auto& hex : path)
        if (isInside(hex)) newPath.push_back(hex);
      std::swap(path, newPath);
    }
    assert(path.size() >= 2);

    // move vertex along path
    for (int j = 0; j < int(path.size()) - 1; ++j) {
      int& curr = array[to1DIndex(path[j].q, path[j].r)];
      int& next = array[to1DIndex(path[j + 1].q, path[j + 1].r)];
      if (next != -1)
        std::swap(points[curr], points[next]);
      else
        points[i] = path[j + 1];
      std::swap(curr, next);
    }
  }

  Eigen::VectorXf toPosition() const {
    // assert(isCorrectState());
    Eigen::VectorXf position(2 * n);
    for (size_t i = 0; i < n; ++i)
      std::tie(position[2 * i], position[2 * i + 1]) = hex2xy(i);
    return position;
  }

  double calcScore(const Problem& problem, bool includeRepulsive = true) const {
    return problem.calcScore(toPosition(), includeRepulsive);
  }

 private:
  double k2;

  // * Only For updateAlongPath
  std::vector<Hex> linedraw(const Hex& a, const Hex& b) const {
    int N = a.distance(b);
    assert(N >= 1);
    if (N == 1) return {a, b};
    int a_s = -a.q - a.r, b_s = -b.q - b.r;
    float a_nudge_q = a.q + 1e-06, a_nudge_r = a.r + 1e-06, a_nudge_s = a_s - 2e-06;
    float b_nudge_q = b.q + 1e-06, b_nudge_r = b.r + 1e-06, b_nudge_s = b_s - 2e-06;
    std::vector<Hex> results;
    float step = 1.0 / std::max(N, 1);
    for (int i = 0; i <= N; ++i)
      results.push_back(Hex::lerp(a_nudge_q, a_nudge_r, a_nudge_s, b_nudge_q, b_nudge_r,
                                  b_nudge_s, step * i));
    return results;
  }

  // * For Util
  inline bool isInside(const Hex& hex) const {
    return 0 <= hex.q && hex.q < int(arraySz) && 0 <= hex.r && hex.r < int(arraySz);
  }

  // * For debugging
  bool isCorrectState() const {
    for (size_t i = 0; i < points.size(); ++i)
      if (array[to1DIndex(points[i].q, points[i].r)] != int(i)) return false;
    return true;
  }

  // * Utility function to convert 2D indices to 1D index
  inline size_t to1DIndex(int q, int r) const { return q * arraySz + r; }
};