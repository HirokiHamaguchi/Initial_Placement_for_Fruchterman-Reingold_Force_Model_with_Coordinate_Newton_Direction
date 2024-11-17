#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "../util/computeDxDy.hpp"
#include "../util/grid.hpp"
#include "../util/hex.hpp"
#include "../util/problem.hpp"
#include "../util/timer.hpp"

// Minimize
// \sum_{e \in adj1} (|\angle (x_{e.u} - x_{e.v})|)
// + \sum_{e \in adj2} (|\angle (x_{e.u} - x_{e.v})|)
// using Simulated Annealing (SA)

struct Circle {
  Circle(size_t n, int seed, std::vector<std::vector<size_t>>& adj)
      : adj1(adj), adj2(n) {
    std::mt19937 g(seed);
    indices.resize(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), g);

    for (size_t i = 0; i < n; ++i) {
      for (size_t j : adj1[i])
        for (size_t k : adj1[j])
          if (i != k && std::find(adj1[i].begin(), adj1[i].end(), k) == adj1[i].end())
            adj2[i].push_back(k);
      std::sort(adj2[i].begin(), adj2[i].end());
      adj2[i].erase(std::unique(adj2[i].begin(), adj2[i].end()), adj2[i].end());
    }
  }

  int angle(int indU, int indV) const {
    int diff = indU - indV;
    if (diff < 0) diff += indices.size();
    assert(0 <= diff && diff < int(indices.size()));
    if (diff > int(indices.size()) / 2) diff = indices.size() - diff;
    assert(0 <= diff && diff <= int(indices.size()) / 2);
    return diff;
  }

  double calcScoreForCircle() {
    double score = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
      for (size_t j : adj1[i]) score += angle(indices[i], indices[j]);
      for (size_t j : adj2[i]) score += angle(indices[i], indices[j]);
    }
    return score / 2;
  }

  // swap indices[u] and indices[v]
  double calcDiff(int u, int v) {
    assert(u != v);
    double diffScore = 0;

    auto updateDiffScore = [&](const std::vector<size_t>& adjList, size_t a, size_t b) {
      for (size_t i : adjList) {
        if (i == b) continue;
        diffScore -= angle(indices[a], indices[i]);
        diffScore += angle(indices[b], indices[i]);
      }
    };

    updateDiffScore(adj1[u], u, v);
    updateDiffScore(adj1[v], v, u);
    updateDiffScore(adj2[u], u, v);
    updateDiffScore(adj2[v], v, u);

    return diffScore;
  }

  std::vector<int> indices;
  std::vector<std::vector<size_t>> adj1;  // theoretical distance = 1
  std::vector<std::vector<size_t>> adj2;  // theoretical distance = 2
};

void solve_circle(const Problem& problem, const int seed,
                  std::vector<Eigen::VectorXf>& positions, Timer& timer) {
  timer.start();

  std::vector<std::vector<size_t>> adj(problem.n);
  for (size_t j = 0; j < problem.m; j++) {
    size_t u = problem.row[j], v = problem.col[j];
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  Circle circle(problem.n, seed, adj);

  const int ITERATIONS = 2 * problem.n * (problem.n * problem.n / problem.m);
  double T0 = +1;
  double T1 = -2;

  double score = circle.calcScoreForCircle();
  double bestScore = score;
  std::vector<int> bestIndices = circle.indices;

  std::mt19937 g(seed);
  std::uniform_int_distribution<int> distVertex(0, problem.n - 1);
  std::uniform_real_distribution<double> distSA(0, 1);

  for (int it = 0; it < ITERATIONS; it++) {
    double T = std::pow(10, T0 + (T1 - T0) * it / ITERATIONS);

    int u = distVertex(g), v = distVertex(g);
    if (u == v) continue;

    double diffScore = circle.calcDiff(u, v);
    if (diffScore <= 0 || distSA(g) < std::exp(-diffScore / T)) {
      std::swap(circle.indices[u], circle.indices[v]);
      score += diffScore;
    }

    if (score < bestScore) {
      bestScore = score;
      bestIndices = circle.indices;
    }
  }

  Eigen::VectorXf finalPos(2 * bestIndices.size());
  for (size_t i = 0; i < bestIndices.size(); ++i) {
    double theta = 2.0 * M_PI * bestIndices[i] / bestIndices.size();
    finalPos[2 * i] = std::cos(theta);
    finalPos[2 * i + 1] = std::sin(theta);
  }

  positions.push_back(finalPos);
  timer.stop();

  double finalScore = circle.calcScoreForCircle();
  circle.indices = bestIndices;
  double bestScore2 = circle.calcScoreForCircle();
  assert(std::abs(finalScore - score) < 1e-5);
  assert(std::abs(bestScore - bestScore2) < 1e-5);
}
