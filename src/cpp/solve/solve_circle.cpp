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
// \sum_{e \in adj1} (||x_{e.u} - x_{e.v}||)
// + \sum_{e \in adj2} (||x_{e.u} - x_{e.v}||)
// using Simulated Annealing (SA)

struct Circle {
  Circle(size_t n, int seed, std::vector<std::vector<int>>& adj) : adj1(adj), adj2(n) {
    std::mt19937 g(seed);
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::shuffle(idx.begin(), idx.end(), g);
    for (size_t i = 0; i < n; ++i) {
      double theta = 2 * M_PI * idx[i] / n;
      positions.emplace_back(std::cos(theta), std::sin(theta));
    }

    for (size_t i = 0; i < n; ++i) {
      for (size_t j : adj1[i])
        for (size_t k : adj1[j])
          if (k != i) adj2[i].push_back(k);
      std::sort(adj2[i].begin(), adj2[i].end());
      adj2[i].erase(std::unique(adj2[i].begin(), adj2[i].end()), adj2[i].end());
    }
  }

  double calcScoreForCircle() {
    double score = 0;
    for (size_t i = 0; i < positions.size(); ++i) {
      for (size_t j : adj1[i])
        score += std::hypot(positions[i].first - positions[j].first,
                            positions[i].second - positions[j].second);
      for (size_t j : adj2[i])
        score += std::hypot(positions[i].first - positions[j].first,
                            positions[i].second - positions[j].second);
    }
    return score / 2;
  }

  // change edge(u, v) to edge(u2, v)
  void calcDiff(int u, int u2, int v, double& diffScore) {
    assert(u != v);
    if (u2 == v) return;
    diffScore -= std::hypot(positions[u].first - positions[u2].first,
                            positions[u].second - positions[u2].second);
    diffScore += std::hypot(positions[v].first - positions[u2].first,
                            positions[v].second - positions[u2].second);
  }

  std::vector<std::pair<double, double>> positions;
  std::vector<std::vector<int>> adj1;  // theoretical distance = 1
  std::vector<std::vector<int>> adj2;  // theoretical distance = 2
};

std::vector<Eigen::VectorXf> solve_circle(const Problem& problem,
                                          const bool measureTime, const int seed,
                                          std::vector<Eigen::VectorXf>& positions,
                                          std::vector<std::pair<double, double>>& hist,
                                          Timer& timer) {
  timer.start();

  std::vector<std::vector<int>> adj(problem.n);
  for (size_t j = 0; j < problem.m; j++) {
    int u = problem.row[j], v = problem.col[j];
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
  Circle circle(problem.n, seed, adj);

  const int ITERATIONS = 0.5 * problem.n * (problem.n * problem.n / problem.m);
  double T0 = +1;
  double T1 = -2;

  double score = circle.calcScoreForCircle();
  double bestScore = score;
  std::vector<std::pair<double, double>> bestPositions = circle.positions;

  std::mt19937 g(seed);
  std::uniform_int_distribution<int> distVertex(0, problem.n - 1);
  std::uniform_real_distribution<double> distSA(0, 1);

  for (int it = 0; it < ITERATIONS; it++) {
    double T = std::pow(10, T0 + (T1 - T0) * it / ITERATIONS);

    int u = distVertex(g), v = distVertex(g);
    if (u == v) continue;

    double diffScore = 0;
    for (auto& u2 : circle.adj1[u]) circle.calcDiff(u, u2, v, diffScore);
    for (auto& v2 : circle.adj1[v]) circle.calcDiff(v, v2, u, diffScore);
    for (auto& u2 : circle.adj2[u]) circle.calcDiff(u, u2, v, diffScore);
    for (auto& v2 : circle.adj2[v]) circle.calcDiff(v, v2, u, diffScore);

    if (diffScore <= 0 || distSA(g) < std::exp(-diffScore / T)) {
      std::swap(circle.positions[u], circle.positions[v]);
      score += diffScore;
    }

    if (score < bestScore) {
      bestScore = score;
      bestPositions = circle.positions;
    }
  }

  Eigen::VectorXf finalPos(2 * bestPositions.size());
  for (size_t i = 0; i < bestPositions.size(); ++i) {
    finalPos[2 * i] = bestPositions[i].first;
    finalPos[2 * i + 1] = bestPositions[i].second;
  }

  positions.push_back(finalPos);
  timer.stop();

  double finalScore = circle.calcScoreForCircle();
  circle.positions = bestPositions;
  double bestScore2 = circle.calcScoreForCircle();
  assert(std::abs(finalScore - score) < 1e-5);
  assert(std::abs(bestScore - bestScore2) < 1e-5);
  assert(!measureTime || positions.size() == 1);

  hist.emplace_back(problem.calcScore(finalPos, true), timer.sec());

  return positions;
}
