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
  Grid(int n, double k, int seed) : n(n), n2(0), k(k) {
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

    std::mt19937 g(seed);
    std::shuffle(points.begin(), points.end(), g);
    points.resize(n);

    array.resize(arraySz * arraySz, -1);
    for (size_t i = 0; i < points.size(); ++i)
      array[to1DIndex(points[i].q, points[i].r)] = i;
  }

  inline std::pair<float, float> hex2xy(double q, double r) const {
    return {q + r / 2.0, r * std::sqrt(3) / 2.0};
  }
  inline std::pair<float, float> hex2xy(size_t i) const {
    return hex2xy(points[i].q, points[i].r);
  }
  Hex xy2hex(float x, float y) {
    float r = y * 2.0 / std::sqrt(3);
    float q = x - r / 2.0;
    return Hex::round(q, r, -q - r);
  }

  void calc_grad_hess(int dq, int dr, float w, float& gx, float& gy, float& hxx,
                      float& hxy, float& hyy) const {
    auto delta = hex2xy(dq, dr);
    float dist = std::hypot(delta.first, delta.second);
    assert(dist > 1e-9);

    // Only use attractive force
    float coeff1 = w * dist / k;
    float coeff2 = w / (dist * k);

    gx += coeff1 * delta.first;
    gy += coeff1 * delta.second;
    hxx += coeff1 + coeff2 * delta.first * delta.first;
    hxy += coeff2 * delta.first * delta.second;
    hyy += coeff1 + coeff2 * delta.second * delta.second;
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

  void swap(int i, const Hex& hexI, const Hex& hexJ) {
    int j = array[to1DIndex(hexJ.q, hexJ.r)];
    if (j == -1) {
      points[i] = hexJ;
      array[to1DIndex(hexI.q, hexI.r)] = -1;
      array[to1DIndex(hexJ.q, hexJ.r)] = i;
    } else {
      std::swap(points[i], points[j]);
      std::swap(array[to1DIndex(hexI.q, hexI.r)], array[to1DIndex(hexJ.q, hexJ.r)]);
    }
  }

  inline bool isInside(const Hex& hex) const {
    return 0 <= hex.q && hex.q < int(arraySz) && 0 <= hex.r && hex.r < int(arraySz);
  }

  // For debugging
  bool isCorrectState() const {
    for (size_t i = 0; i < points.size(); ++i)
      if (array[to1DIndex(points[i].q, points[i].r)] != int(i)) return false;
    int cnt = std::count_if(array.begin(), array.end(), [](int x) { return x != -1; });
    return cnt == int(points.size());
  }

 private:
  // Utility function to convert 2D indices to 1D index
  inline size_t to1DIndex(int q, int r) const { return q * arraySz + r; }
};