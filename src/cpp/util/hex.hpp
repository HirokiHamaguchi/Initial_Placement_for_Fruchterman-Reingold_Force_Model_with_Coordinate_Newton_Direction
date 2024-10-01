#pragma once

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "problem.hpp"

/**
 * @brief Hex class
 * @ref https://www.redblobgames.com/grids/hexagons/
 */
class Hex {
 public:
  int q, r;

  Hex() : q(-INT_MAX), r(-INT_MAX) {}
  Hex(int q, int r) : q(q), r(r) {}

  Hex operator+(const Hex& other) const { return Hex(q + other.q, r + other.r); }

  Hex operator-(const Hex& other) const { return Hex(q - other.q, r - other.r); }

  bool operator==(const Hex& other) const { return q == other.q && r == other.r; }

  bool operator!=(const Hex& other) const { return !(*this == other); }

  int length() const { return (std::abs(q) + std::abs(r) + std::abs(-q - r)) / 2; }

  int distance(const Hex& other) const { return (*this - other).length(); }

  friend std::ostream& operator<<(std::ostream& os, const Hex& hex) {
    os << "Hex(" << hex.q << ", " << hex.r << ")";
    return os;
  }

  static Hex round(float q, float r, float s) {
    int qi = std::round(q);
    int ri = std::round(r);
    int si = std::round(s);
    float q_diff = std::abs(qi - q);
    float r_diff = std::abs(ri - r);
    float s_diff = std::abs(si - s);
    if (q_diff > r_diff && q_diff > s_diff) {
      qi = -ri - si;
    } else if (r_diff > s_diff) {
      ri = -qi - si;
    }
    return Hex(qi, ri);
  }

  static Hex lerp(float a_q, float a_r, float a_s, float b_q, float b_r, float b_s,
                  float t) {
    return Hex::round(a_q * (1.0 - t) + b_q * t, a_r * (1.0 - t) + b_r * t,
                      a_s * (1.0 - t) + b_s * t);
  }

  static bool equal_hex_array(const std::vector<Hex>& a, const std::vector<Hex>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i] != b[i]) return false;
    return true;
  }
};

struct Grid {
  int n;   // number of vertices
  int n2;  // length of the side of the hexagon
  double k;
  std::vector<Hex> points;
  std::vector<std::vector<int>> array;

 private:
  double k2;
  double kForHex;

 public:
  Grid(int n, double k) : n(n), n2(0), k(k), k2(std::pow(k, 2)), kForHex(k) {
    int hexSize = 2 * n;
    while (3 * n2 * n2 + 3 * n2 + 1 < hexSize) n2++;
    for (int r = 0; r <= 2 * n2; ++r) {
      for (int q = 0; q <= 2 * n2; ++q) {
        if (r + q < n2 || r + q > 3 * n2) continue;
        points.emplace_back(q, r);
      }
    }
    assert(int(points.size()) >= hexSize);

    std::mt19937 g(0);
    std::shuffle(points.begin(), points.end(), g);
    points.resize(n);

    array.resize(2 * n2 + 1, std::vector<int>(2 * n2 + 1, -1));
    for (size_t i = 0; i < points.size(); ++i) array[points[i].q][points[i].r] = i;
  }

  void setKForHex(double kForHex) { this->kForHex = kForHex; }

  inline bool isInside(const Hex& hex) const {
    return 0 <= hex.q && hex.q < 2 * n2 + 1 && 0 <= hex.r && hex.r < 2 * n2 + 1;
  }
  inline std::pair<float, float> hex2xy(double q, double r) const {
    return {kForHex * (q + r / 2.0), kForHex * (r * std::sqrt(3) / 2.0)};
  }
  inline std::pair<float, float> hex2xy(int i) const {
    return hex2xy(points[i].q, points[i].r);
  }
  Hex xy2hex(float x, float y) {
    float r = y * 2.0 / (kForHex * std::sqrt(3));
    float q = x / kForHex - r / 2.0;
    return Hex::round(q, r, -q - r);
  }

  bool isCorrectState() const {
    for (size_t i = 0; i < points.size(); ++i)
      if (array[points[i].q][points[i].r] != int(i)) return false;
    return true;
  }

  void calc_grad_hess(int dq, int dr, double w, double& gx, double& gy, double& hxx,
                      double& hxy, double& hyy) const {
    auto delta = hex2xy(dq, dr);
    double dist = std::hypot(delta.first, delta.second);
    assert(dist > 1e-9);

    // * Method 1: only use attractive force
    double coeff1 = w * dist / k;
    double coeff2 = w / (dist * k);

    // * Method 2: use both attractive and repulsive forces
    // double d2 = std::pow(dist, 2);
    // double d4 = std::pow(d2, 2);
    // double coeff1 = w * dist / k - k2 / d2;
    // double coeff2 = w / (dist * k) + 2 * k2 / d4;

    gx += coeff1 * delta.first;
    gy += coeff1 * delta.second;
    hxx += coeff1 + coeff2 * delta.first * delta.first;
    hxy += coeff2 * delta.first * delta.second;
    hyy += coeff1 + coeff2 * delta.second * delta.second;
  }

  std::vector<Hex> linedraw(const Hex& a, const Hex& b) const {
    int N = a.distance(b);
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

  double calcScore(const Problem& problem, bool includeRepulsive = true) const {
    double score = 0.0;

    if (includeRepulsive) {
      double k2 = std::pow(k, 2);
      for (size_t u = 0; u < problem.n; ++u) {
        for (size_t v = u + 1; v < problem.n; ++v) {
          Hex hex = points[u] - points[v];
          auto [dx, dy] = hex2xy(hex.q, hex.r);
          double d = std::hypot(dx, dy);
          assert(d > 1e-9);
          score -= k2 * std::log(d);
        }
      }
    }

    for (size_t i = 0; i < problem.m; ++i) {
      size_t u = problem.row[i];
      size_t v = problem.col[i];
      assert(u < v);
      double w = problem.data[i];
      Hex hex = points[u] - points[v];
      auto [dx, dy] = hex2xy(hex.q, hex.r);
      double d = std::hypot(dx, dy);
      score += w * std::pow(d, 3) / (3.0 * problem.k);
    }

    return score;
  }

  double computeKForHex(const Problem& problem) {
    // Minimize_x x^3 score_a - k^2 n(n-1) \log(x)
    // where score_a = \sum_{i < j} w_{ij} d_{ij}^3 / (3k)
    double score_a = calcScore(problem, false);
    double coeff_r = std::pow(k, 2) * n * (n - 1);

    // Minimize f(x) = x^3 score_a - coeff_r \log(x) : convex
    // f'(x) = 3x^2 score_a - coeff_r / x
    // f''(x) = 6x score_a + coeff_r / x^2
    double x = 1.0;
    for (int iter = 0; iter < 10; ++iter) {
      double df = 3 * std::pow(x, 2) * score_a - coeff_r / x;
      double ddf = 6 * x * score_a + coeff_r / std::pow(x, 2);
      double dx = df / ddf;
      if (std::abs(dx) < 1e-9) break;
      x -= dx;
    }
    dbg(x);
    setKForHex(kForHex * x);
    return kForHex * x;
  }

  Eigen::VectorXf toPosition() const {
    assert(isCorrectState());
    Eigen::VectorXf position(2 * n);
    for (int i = 0; i < n; ++i)
      std::tie(position[2 * i], position[2 * i + 1]) = hex2xy(i);
    return position;
  }
};
