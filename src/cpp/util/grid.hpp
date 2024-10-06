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

template <class E>
struct csr {
  vector<int> start;
  vector<E> elist;
  explicit csr(int n, const vector<pair<int, E>>& edges)
      : start(n + 1), elist(edges.size()) {
    for (auto e : edges) start[e.first + 1]++;
    for (int i = 1; i <= n; i++) start[i] += start[i - 1];
    auto counter = start;
    for (auto e : edges) elist[counter[e.first]++] = e.second;
  }
};
template <class T>
struct simple_queue {
  vector<T> payload;
  int pos = 0;
  void reserve(int n) { payload.reserve(n); }
  int size() const { return int(payload.size()) - pos; }
  bool empty() const { return pos == int(payload.size()); }
  void push(const T& t) { payload.push_back(t); }
  T& front() { return payload[pos]; }
  void clear() {
    payload.clear();
    pos = 0;
  }
  void pop() { pos++; }
};
template <class Cap = long long, class Cost = long long>
struct mcf_graph {
  mcf_graph() {};
  explicit mcf_graph(int n) : _n(n) {};
  struct edge {
    int from, to;
    Cap cap, flow;
    Cost cost;
  };
  int add_edge(int from, int to, Cap cap, Cost cost) {
    assert(0 <= from && from < _n);
    assert(0 <= to && to < _n);
    assert(0 <= cap);
    assert(0 <= cost);
    int m = int(_edges.size());
    _edges.push_back({from, to, cap, 0, cost});
    return m;
  }
  edge get_edge(int i) {
    int m = int(_edges.size());
    assert(0 <= i && i < m);
    return _edges[i];
  }
  vector<edge> edges() { return _edges; }
  pair<Cap, Cost> flow(int s, int t) { return flow(s, t, numeric_limits<Cap>::max()); }
  pair<Cap, Cost> flow(int s, int t, Cap flow_limit) {
    return slope(s, t, flow_limit).back();
  }
  // 返り値の最後の要素は y = min(x, flow_limit)として (y,g(y))
  vector<pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) {
    assert(0 <= s && s < _n);
    assert(0 <= t && t < _n);
    assert(s != t);
    int m = int(_edges.size());
    vector<int> edge_idx(m);
    auto g = [&]() {
      vector<int> degree(_n), redge_idx(m);
      vector<pair<int, _edge>> elist;
      elist.reserve(2 * m);
      for (int i = 0; i < m; i++) {
        auto e = _edges[i];
        edge_idx[i] = degree[e.from]++;
        redge_idx[i] = degree[e.to]++;
        elist.push_back({e.from, {e.to, -1, e.cap - e.flow, e.cost}});
        elist.push_back({e.to, {e.from, -1, e.flow, -e.cost}});
      }
      auto _g = csr<_edge>(_n, elist);
      for (int i = 0; i < m; i++) {
        auto e = _edges[i];
        edge_idx[i] += _g.start[e.from];
        redge_idx[i] += _g.start[e.to];
        _g.elist[edge_idx[i]].rev = redge_idx[i];
        _g.elist[redge_idx[i]].rev = edge_idx[i];
      }
      return _g;
    }();
    auto result = slope(g, s, t, flow_limit);
    for (int i = 0; i < m; i++) {
      auto e = g.elist[edge_idx[i]];
      _edges[i].flow = _edges[i].cap - e.cap;
    }
    return result;
  }

 private:
  int _n;
  vector<edge> _edges;
  struct _edge {
    int to, rev;
    Cap cap;
    Cost cost;
  };
  vector<pair<Cap, Cost>> slope(csr<_edge>& g, int s, int t, Cap flow_limit) {
    // variants (C = maxcost): -(n-1)C <= dual[s] <= dual[i] <= dual[t] = 0 reduced cost
    // (= e.cost + dual[e.from] - dual[e.to]) >= 0 for all edge
    vector<pair<Cost, Cost>> dual_dist(_n);  // dual_dist[i] = (dual[i], dist[i])
    vector<int> prev_e(_n);
    vector<bool> vis(_n);
    struct Q {
      Cost key;
      int to;
      bool operator<(Q r) const { return key > r.key; }
    };
    vector<int> que_min;
    vector<Q> que;
    auto dual_ref = [&]() {
      for (int i = 0; i < _n; i++) {
        dual_dist[i].second = numeric_limits<Cost>::max();
      }
      fill(vis.begin(), vis.end(), false);
      que_min.clear();
      que.clear(); /* que[0..heap_r) was heapified */
      size_t heap_r = 0;
      dual_dist[s].second = 0;
      que_min.push_back(s);
      while (!que_min.empty() || !que.empty()) {
        int v;
        if (!que_min.empty()) {
          v = que_min.back();
          que_min.pop_back();
        } else {
          while (heap_r < que.size()) {
            heap_r++;
            push_heap(que.begin(), que.begin() + heap_r);
          }
          v = que.front().to;
          pop_heap(que.begin(), que.end());
          que.pop_back();
          heap_r--;
        }
        if (vis[v]) continue;
        vis[v] = true;
        if (v == t) break;
        // dist[v] = shortest(s, v) + dual[s] - dual[v],  dist[v] >= 0 (all reduced cost
        // are positive), dist[v] <= (n-1)C
        Cost dual_v = dual_dist[v].first, dist_v = dual_dist[v].second;
        for (int i = g.start[v]; i < g.start[v + 1]; i++) {
          auto e = g.elist[i];
          if (!e.cap)
            continue;  // |-dual[e.to] + dual[v]| <= (n-1)C  cost <= C - -(n-1)C + 0 =
                       // nC
          Cost cost = e.cost - dual_dist[e.to].first + dual_v;
          if (dual_dist[e.to].second - dist_v > cost) {
            Cost dist_to = dist_v + cost;
            dual_dist[e.to].second = dist_to;
            prev_e[e.to] = e.rev;
            if (dist_to == dist_v)
              que_min.push_back(e.to);
            else
              que.push_back(Q{dist_to, e.to});
          }
        }
      }
      if (!vis[t]) return false;
      for (int v = 0; v < _n; v++) {
        if (!vis[v]) continue;
        dual_dist[v].first -= dual_dist[t].second - dual_dist[v].second;
      }
      return true;
    };
    Cap flow = 0;
    Cost cost = 0, prev_cost_per_flow = -1;
    vector<pair<Cap, Cost>> result = {{Cap(0), Cost(0)}};
    while (flow < flow_limit) {
      if (!dual_ref()) break;
      Cap c = flow_limit - flow;
      for (int v = t; v != s; v = g.elist[prev_e[v]].to)
        c = min(c, g.elist[g.elist[prev_e[v]].rev].cap);
      for (int v = t; v != s; v = g.elist[prev_e[v]].to) {
        auto& e = g.elist[prev_e[v]];
        e.cap += c;
        g.elist[e.rev].cap -= c;
      }
      Cost d = -dual_dist[s].first;
      flow += c;
      cost += c * d;
      if (prev_cost_per_flow == d) result.pop_back();
      result.push_back({flow, cost});
      prev_cost_per_flow = d;
    }
    return result;
  }
};

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

    std::vector<Hex> newPoints(n);
    std::vector<int> indexes(n);
    std::iota(indexes.begin(), indexes.end(), 0);
    std::shuffle(indexes.begin(), indexes.end(), std::mt19937(loopCnt));
    for (int i : indexes) {
      double x = newPos[2 * i], y = newPos[2 * i + 1];
      Hex new_v = xy2hex(x, y);
      int cnt = 0;
      for (const auto& deltaHex : deltaHexList) {
        cnt++;
        Hex hex = new_v + deltaHex;
        if (!isInside(hex)) continue;
        if (array[hex.q][hex.r] == -1) {
          newPoints[i] = hex;
          array[hex.q][hex.r] = i;
          break;
        }
      }
      assert(newPoints[i].q != -INT_MAX);
    }
    std::swap(points, newPoints);
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
