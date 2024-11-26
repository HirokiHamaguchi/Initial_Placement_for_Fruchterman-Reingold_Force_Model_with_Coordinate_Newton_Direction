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

  static bool equal_hex_array(const std::vector<Hex>& a, const std::vector<Hex>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a[i] != b[i]) return false;
    return true;
  }
};
