#pragma once

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

/**
 * @brief Hex class
 * @ref https://www.redblobgames.com/grids/hexagons/
 */
class Hex {
 public:
  int q, r, s;

  Hex(int q, int r, int s) : q(q), r(r), s(s) { assert(q + r + s == 0); }

  Hex operator+(const Hex& other) const {
    return Hex(q + other.q, r + other.r, s + other.s);
  }

  Hex operator-(const Hex& other) const {
    return Hex(q - other.q, r - other.r, s - other.s);
  }

  bool operator==(const Hex& other) const {
    return q == other.q && r == other.r && s == other.s;
  }

  bool operator!=(const Hex& other) const { return !(*this == other); }

  std::tuple<float, float> hex2xy(float k) const {
    return std::make_tuple(k * (q + r / 2.0), k * (r * std::sqrt(3) / 2.0));
  }

  int length() const { return (std::abs(q) + std::abs(r) + std::abs(s)) / 2; }

  int distance(const Hex& other) const { return (*this - other).length(); }

  friend std::ostream& operator<<(std::ostream& os, const Hex& hex) {
    os << "Hex(" << hex.q << ", " << hex.r << ", " << hex.s << ")";
    return os;
  }

  static Hex xy2hex(float x, float y, float k) {
    float r = std::round(y * 2.0 / (k * std::sqrt(3)));
    float q = std::round(x / k - r / 2.0);
    return Hex::round(q, r, -q - r);
  }

  static Hex round(float q, float r, float s) {
    int qi = std::round(q);
    int ri = std::round(r);
    int si = std::round(s);
    float q_diff = std::abs(qi - q);
    float r_diff = std::abs(ri - r);
    float s_diff = std::abs(si - s);
    if (q_diff > r_diff && q_diff > s_diff)
      qi = -ri - si;
    else if (r_diff > s_diff)
      ri = -qi - si;
    else
      si = -qi - ri;
    return Hex(qi, ri, si);
  }

  static Hex lerp(float a_q, float a_r, float a_s, float b_q, float b_r, float b_s,
                  float t) {
    return Hex::round(a_q * (1.0 - t) + b_q * t, a_r * (1.0 - t) + b_r * t,
                      a_s * (1.0 - t) + b_s * t);
  }

  static std::vector<Hex> linedraw(const Hex& a, const Hex& b) {
    int N = a.distance(b);
    float a_nudge_q = a.q + 1e-06, a_nudge_r = a.r + 1e-06, a_nudge_s = a.s - 2e-06;
    float b_nudge_q = b.q + 1e-06, b_nudge_r = b.r + 1e-06, b_nudge_s = b.s - 2e-06;
    std::vector<Hex> results;
    float step = 1.0 / std::max(N, 1);
    for (int i = 0; i <= N; ++i) {
      results.push_back(Hex::lerp(a_nudge_q, a_nudge_r, a_nudge_s, b_nudge_q, b_nudge_r,
                                  b_nudge_s, step * i));
    }
    return results;
  }

  static bool equal_hex_array(const std::vector<Hex>& a, const std::vector<Hex>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i] != b[i]) return false;
    return true;
  }
};
