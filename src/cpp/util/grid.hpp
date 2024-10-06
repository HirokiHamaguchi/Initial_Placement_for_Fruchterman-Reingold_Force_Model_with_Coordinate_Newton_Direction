#pragma once

#include <algorithm>
#include <atcoder/mincostflow>
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
  int n;   // number of vertices
  int n2;  // length of the side of the hexagon
  double k;
  std::vector<Hex> points;
  std::vector<std::vector<int>> array;

 public:
  Grid(int n, double k) : n(n), n2(0), k(k), k2(std::pow(k, 2)) {
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

    initializeDeltaHexList();
  }

  inline std::pair<float, float> hex2xy(double q, double r) const {
    return {k * (q + r / 2.0), k * (r * std::sqrt(3) / 2.0)};
  }
  inline std::pair<float, float> hex2xy(int i) const {
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

  void updateAlongPath(int i, const Hex& new_v) {
    std::vector<Hex> path;
    for (auto& hex : linedraw(points[i], new_v))
      if (isInside(hex)) path.push_back(hex);
    assert(path.size() >= 2);

    // move vertex along path
    for (int j = 0; j < int(path.size()) - 1; ++j) {
      int& curr = array[path[j].q][path[j].r];
      int& next = array[path[j + 1].q][path[j + 1].r];
      if (next != -1)
        std::swap(points[curr], points[next]);
      else
        points[i] = path[j + 1];
      std::swap(curr, next);
    }
  }

  bool updateToNewPos(const Eigen::VectorXf& _newPos, const int loopCnt) {
    for (int i = 0; i < n; ++i) array[points[i].q][points[i].r] = -1;

    Eigen::VectorXf newPos = _newPos;
    float cx = 0.0, cy = 0.0;
    for (int i = 0; i < n; ++i) {
      cx += newPos[2 * i], cy += newPos[2 * i + 1];
    }
    cx /= n, cy /= n;
    auto [originX, originY] = hex2xy(n2, n2);
    for (int i = 0; i < n; ++i) {
      newPos[2 * i] += originX - cx;
      newPos[2 * i + 1] += originY - cy;
    }

    atcoder::mcf_graph<long long, long long> graph(n + (2 * n2 + 1) * (2 * n2 + 1) + 2);
    int s = n + (2 * n2 + 1) * (2 * n2 + 1), t = s + 1;
    // source -> vertex
    for (int i = 0; i < n; ++i) graph.add_edge(s, i, 1, 0);
    // vertex -> hex
    for (int i = 0; i < n; ++i) {
      double x = newPos[2 * i], y = newPos[2 * i + 1];
      Hex hex = xy2hex(x, y);
      for (int dq = -50; dq <= 50; dq++) {
        for (int dr = -50; dr <= 50; dr++) {
          if (std::abs(dq) + std::abs(dr) > 50) continue;
          Hex newHex = hex + Hex(dq, dr);
          if (!isInside(newHex)) continue;
          auto [x2, y2] = hex2xy(newHex.q, newHex.r);
          int hexIdx = n + (newHex.q * (2 * n2 + 1) + newHex.r);
          long long cost =
              std::llround(std::pow(x - x2, 2) + std::pow(y - y2, 2) * 1000000);
          graph.add_edge(i, hexIdx, 1, cost);
        }
      }
    }
    // hex -> sink
    for (int q = 0; q <= 2 * n2; ++q) {
      for (int r = 0; r <= 2 * n2; ++r) {
        graph.add_edge(n + (q * (2 * n2 + 1) + r), t, 1, 0);
      }
    }

    auto [flow, cost] = graph.flow(s, t, n);
    dbg(flow, cost);
    if (flow != n) assert(false);

    std::vector<Hex> newPoints(n);
    auto edges = graph.edges();
    for (auto& e : edges) {
      if (e.from < n && e.flow == 1) {
        int i = e.from;
        int hexIdx = e.to - n;
        assert(newPoints[i].q == -INT_MAX && newPoints[i].r == -INT_MAX);
        newPoints[i] = Hex(hexIdx / (2 * n2 + 1), hexIdx % (2 * n2 + 1));
      }
    }

    assert(std::all_of(newPoints.begin(), newPoints.end(), [&](const Hex& hex) {
      return hex.q != -INT_MAX && hex.r != -INT_MAX;
    }));
    std::swap(points, newPoints);
    for (int i = 0; i < n; ++i) array[points[i].q][points[i].r] = i;
    return points != newPoints;
  }

  Eigen::VectorXf toPosition() const {
    assert(isCorrectState());
    Eigen::VectorXf position(2 * n);
    for (int i = 0; i < n; ++i)
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
    return 0 <= hex.q && hex.q < 2 * n2 + 1 && 0 <= hex.r && hex.r < 2 * n2 + 1;
  }

  // * For updateToNewPos
  std::vector<Hex> deltaHexList;
  void initializeDeltaHexList() {
    assert(deltaHexList.empty());
    int maxDist = 2 * n2 + 1;
    for (int dq = -maxDist; dq <= +maxDist; ++dq)
      for (int dr = -maxDist; dr <= +maxDist; ++dr) deltaHexList.emplace_back(dq, dr);
    std::sort(deltaHexList.begin(), deltaHexList.end(),
              [&](const Hex& a, const Hex& b) {
                auto [x1, y1] = hex2xy(a.q, a.r);
                auto [x2, y2] = hex2xy(b.q, b.r);
                return (x1 * x1 + y1 * y1) < (x2 * x2 + y2 * y2);
              });
  }

  // * For debugging
  bool isCorrectState() const {
    for (size_t i = 0; i < points.size(); ++i)
      if (array[points[i].q][points[i].r] != int(i)) return false;
    return true;
  }
};
