#include "include/LBFGS.h"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"
#include "util/timer.hpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  RS_FR,      // Init:RS / Optimize:FR
  RS_L_BFGS,  // Init:RS / Optimize:L_BFGS
};

std::string MethodStr[4] = {"FR", "L_BFGS", "RS_FR", "RS_L_BFGS"};

std::pair<std::vector<std::pair<double, double>>, std::vector<Eigen::VectorXf>> solve(
    const Method method, const Problem& problem, const bool measureTime,
    const int seed) {
  std::vector<std::pair<double, double>> hist;
  std::vector<Eigen::VectorXf> positions;

  Timer timer;

  if (method == RS_FR || method == RS_L_BFGS) {
    solve_init(problem, measureTime, seed, positions, hist, timer);
  } else {
    timer.start();
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) position[i] = std::abs(position[i]);
    positions.push_back(position);
    timer.stop();
  }

  std::vector<Eigen::VectorXf> positions2;
  if (method == L_BFGS || method == RS_L_BFGS) {
    solve_LBFGS<FunctionFR>(problem, positions, hist, timer);
  } else if (method == FR || method == RS_FR) {
    solve_FR(problem, positions, hist, timer);
  }

  assert(!timer.isMeasuring);

  return {hist, positions};
}
