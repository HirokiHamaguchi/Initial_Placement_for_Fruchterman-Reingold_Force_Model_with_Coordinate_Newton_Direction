#include "../include/LBFGS.h"
#include "../util/function.hpp"
#include "../util/problem.hpp"
#include "../util/timer.hpp"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  SN_FR,      // Init:SN / Optimize:FR
  SN_L_BFGS,  // Init:SN / Optimize:L_BFGS
};

std::string MethodStr[4] = {"FR", "L_BFGS", "SN_FR", "SN_L_BFGS"};

std::pair<std::vector<std::pair<double, double>>, std::vector<Eigen::VectorXf>> solve(
    const Method method, const Problem& problem, const bool measureTime, const int seed,
    const int MAX_ITER) {
  std::vector<std::pair<double, double>> hist;
  std::vector<Eigen::VectorXf> positions;

  Timer timer;

  if (method == SN_FR || method == SN_L_BFGS) {
    solve_init(problem, measureTime, seed, positions, hist, timer);
  } else {
    timer.start();
    std::srand(seed);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < int(position.size()); i++) position[i] = std::abs(position[i]);
    positions.push_back(position);
    timer.stop();
  }

  std::vector<Eigen::VectorXf> positions2;
  if (method == L_BFGS || method == SN_L_BFGS) {
    solve_LBFGS<FunctionFR>(problem, positions, hist, timer, MAX_ITER);
  } else if (method == FR || method == SN_FR) {
    solve_FR(problem, positions, hist, timer, MAX_ITER);
  }
  assert(!timer.isMeasuring);

  return {hist, positions};
}
