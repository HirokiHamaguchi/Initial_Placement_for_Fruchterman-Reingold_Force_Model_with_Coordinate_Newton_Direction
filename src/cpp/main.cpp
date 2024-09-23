#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <vector>

#include "util/problem.hpp"

std::vector<std::pair<double, double>> toHexagon(
    const std::vector<std::pair<double, double>>& _pos, double k) {
  std::vector<std::pair<double, double>> pos = _pos;
  int n = pos.size();
  int sqrtN = static_cast<int>(std::sqrt(n));
  double invSqrt3 = 1 / std::sqrt(3.0);

  std::vector<int> args(n);
  std::iota(args.begin(), args.end(), 0);
  std::sort(args.begin(), args.end(),
            [&pos](int a, int b) { return pos[a].first < pos[b].first; });
  for (int i = 0; i < n; i++) pos[args[i]].first = i / sqrtN;

  std::sort(args.begin(), args.end(), [&pos](int a, int b) { return pos[a] < pos[b]; });

  for (int i = 0; i < n; ++i) {
    pos[args[i]].second =
        (i % sqrtN) * k + (static_cast<int>(pos[args[i]].first) % 2) * 0.5 * k;
    pos[args[i]].first *= 1.5 * invSqrt3 * k;
  }

  return pos;
}

std::vector<std::vector<std::pair<double, double>>> solve(
    const Problem& problem, std::vector<std::pair<double, double>>& pos) {
  int n = problem.n;
  double k = problem.k;

  pos = toHexagon(pos, k);
  std::vector<std::vector<std::pair<double, double>>> positions = {pos};

  std::vector<double> distance(n);

  auto t0 = std::chrono::high_resolution_clock::now();

  const int iterations = 1000;
  for (int iter = 0; iter < iterations; ++iter) {
    std::vector<std::pair<double, double>> displacement(n, {0.0, 0.0});

    for (size_t i = 0; i < problem.m; i++) {
      int u = problem.row[i];
      int v = problem.col[i];
      double delta_x = pos[u].first - pos[v].first;
      double delta_y = pos[u].second - pos[v].second;
      double dist = std::hypot(delta_x, delta_y);
      double factor = -problem.data[i] * dist / k;
      displacement[u].first += factor * delta_x;
      displacement[u].second += factor * delta_y;
      displacement[v].first -= factor * delta_x;
      displacement[v].second -= factor * delta_y;
    }

    for (int i = 0; i < n; ++i) {
      double norm = std::hypot(displacement[i].first, displacement[i].second);
      if (norm > 0) {
        displacement[i].first /= norm;
        displacement[i].second /= norm;
      }
      displacement[i].first *= 2 * k;
      displacement[i].second *= 2 * k;
    }

    for (int i = 0; i < n; ++i) {
      pos[i].first += displacement[i].first;
      pos[i].second += displacement[i].second;
    }

    pos = toHexagon(pos, k);
    // dbg(problem.calcScore(pos, false));
    positions.push_back(pos);
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = t1 - t0;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

  return positions;
}

int main() {
  // Step1 : Read input
  Problem problem("494_bus");

  // Step2 : decide initial position
  std::vector<std::pair<double, double>> position(problem.n);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(0, 1);
  for (int i = 0; i < int(problem.n); ++i) {
    position[i] = {dist(gen), dist(gen)};
  }
  
  // Step3 : Solve
  std::vector<std::vector<std::pair<double, double>>> positions =
      solve(problem, position);

  // Step4 : Output
  problem.printOutput(positions);

  return 0;
}

// std::pair<double, double> axialCoord2Euclidean(int q, int r) {
//   return {q + r / 2.0, r * std::sqrt(3) / 2.0};
// }

// std::vector<std::pair<double, double>> hexagon(int n, double k) {
//   std::vector<std::pair<double, double>> ret = {{0, 0}};
//   if (n == 1) return ret;

//   std::vector<std::tuple<int, int, int>> diffs = {{0, 1, -1}, {-1, 1, 0}, {-1, 0, 1},
//                                                   {0, -1, 1}, {1, -1, 0}, {1, 0,
//                                                   -1}};

//   int distFromCenter = 0;
//   while (true) {
//     int q = distFromCenter, r = -distFromCenter, s = 0;
//     for (int i = 0; i < 6; ++i) {
//       for (int j = 0; j < distFromCenter; ++j) {
//         ret.push_back(axialCoord2Euclidean(k * q, k * r));
//         if (int(ret.size()) == n) return ret;
//         q += std::get<0>(diffs[i]);
//         r += std::get<1>(diffs[i]);
//         s += std::get<2>(diffs[i]);
//       }
//     }
//     distFromCenter += 1;
//   }
// }
