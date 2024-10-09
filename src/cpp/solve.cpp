#include "include/LBFGS.h"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  RS_FR,      // Init:RS / Optimize:FR
  RS_L_BFGS,  // Init:RS / Optimize:L_BFGS
};

std::string MethodStr[4] = {"FR", "L_BFGS", "RS_FR", "RS_L_BFGS"};

std::tuple<std::vector<double>, std::vector<Eigen::VectorXf>, double> solve(
    const Method method, const Problem& problem, const bool measureTime,
    const int seed) {
  std::vector<double> hist;
  std::vector<Eigen::VectorXf> positions;

  auto t0 = std::chrono::high_resolution_clock::now();
  if (method == RS_FR || method == RS_L_BFGS) {
    positions = solve_init(problem, measureTime, seed);
  } else {
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) position[i] = std::abs(position[i]);
    positions.push_back(position);
  }
  std::vector<Eigen::VectorXf> positions2;
  if (method == L_BFGS || method == RS_L_BFGS) {
    std::tie(hist, positions2) = solve_LBFGS<FunctionFR>(problem, positions.back());
  } else if (method == FR || method == RS_FR) {
    positions2 = solve_FR(problem, positions.back());
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  double elapsedTime = std::chrono::duration<double>(t1 - t0).count();

  if (method == FR || method == RS_FR)
    for (auto& pos : positions2) hist.push_back(problem.calcScore(pos));
  positions.insert(positions.end(), positions2.begin(), positions2.end());
  return {hist, positions, elapsedTime};
}
