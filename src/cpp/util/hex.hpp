#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <vector>

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

  std::pair<float, float> hex2xy(float k) const {
    return {k * (q + r / 2.0), k * (r * std::sqrt(3) / 2.0)};
  }

  int length() const { return (std::abs(q) + std::abs(r) + std::abs(-q - r)) / 2; }

  int distance(const Hex& other) const { return (*this - other).length(); }

  friend std::ostream& operator<<(std::ostream& os, const Hex& hex) {
    os << "Hex(" << hex.q << ", " << hex.r << ")";
    return os;
  }

  static Hex xy2hex(float x, float y, float k) {
    float r = y * 2.0 / (k * std::sqrt(3));
    float q = x / k - r / 2.0;
    return Hex::round(q, r, -q - r);
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

class Grid {
 public:
  int n;   // number of vertices
  int n2;  // length of the side of the hexagon
  float k;
  std::vector<Hex> points;
  std::vector<std::vector<int>> array;
  std::vector<std::tuple<double, double, double, double, double>> distStore;
  // std::map<std::pair<int, int>, float> logDistStore;

  Grid(int n, float k) : n(n), n2(0), k(k) {
    while (3 * n2 * n2 + 3 * n2 + 1 < n) n2++;

    for (int r = 0; r <= 2 * n2; ++r) {
      for (int q = 0; q <= 2 * n2; ++q) {
        if (r + q < n2 || r + q > 3 * n2) continue;
        points.emplace_back(q, r);
      }
    }
    assert(int(points.size()) >= n);

    std::mt19937 g(0);
    std::shuffle(points.begin(), points.end(), g);
    points.resize(n);

    array.resize(2 * n2 + 1, std::vector<int>(2 * n2 + 1, -1));
    for (size_t i = 0; i < points.size(); ++i) array[points[i].q][points[i].r] = i;

    distStore.resize(400);
    for (int q = -10; q < 10; ++q) {
      for (int r = -10; r < 10; ++r) {
        distStore[20 * (q + 10) + (r + 10)] = calc_grad_hess<true>(q, r, 1.0);
        // logDistStore[20 * (q + 10) + (r + 10)] =
        //     (q == 0 && r == 0) ? std::numeric_limits<float>::infinity()
        //                        : std::log(distStore[{q, r}].second);
      }
    }
  }

  bool isCorrectState() const {
    for (size_t i = 0; i < points.size(); ++i)
      if (array[points[i].q][points[i].r] != int(i)) return false;
    return true;
  }

  template <bool init = false>
  std::tuple<double, double, double, double, double> calc_grad_hess(int dq, int dr,
                                                                    double w) const {
    if (!init && abs(dq) < 10 && abs(dr) < 10) {
      auto [gx, gy, hxx, hxy, hyy] = distStore[20 * (dq + 10) + (dr + 10)];
      return {w * gx, w * gy, w * hxx, w * hxy, w * hyy};
    } else {
      auto delta = Hex(dq, dr).hex2xy(k);
      double dist = std::hypot(delta.first, delta.second);
      double coeff1 = dist / k;
      double coeff2 = 1 / (dist * k);
      double gx = coeff1 * delta.first;
      double gy = coeff1 * delta.second;
      double hxx = coeff1 + coeff2 * delta.first * delta.first;
      double hxy = coeff2 * delta.first * delta.second;
      double hyy = coeff1 + coeff2 * delta.second * delta.second;
      return {w * gx, w * gy, w * hxx, w * hxy, w * hyy};
    }
  }

  bool isInside(const Hex& hex) const {
    return 0 <= hex.q && hex.q < 2 * n2 + 1 && 0 <= hex.r && hex.r < 2 * n2 + 1;
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
    assert(std::all_of(results.begin(), results.end(),
                       [this](const Hex& hex) { return isInside(hex); }));
    return results;
  }

  double calcScore(const Problem& problem, bool includeRepulsive = true) const {
    double score = 0.0;

    if (includeRepulsive) {
      double k2 = std::pow(k, 2);
      for (size_t u = 0; u < problem.n; ++u) {
        for (size_t v = u + 1; v < problem.n; ++v) {
          auto [dx, dy] = (points[u] - points[v]).hex2xy(k);
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
      auto [dx, dy] = (points[u] - points[v]).hex2xy(k);
      double d = std::hypot(dx, dy);
      score += w * std::pow(d, 3) / (3.0 * problem.k);
    }

    return score;
  }
};
