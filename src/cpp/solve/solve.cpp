#include "../include/LBFGS.h"
#include "../util/function.hpp"
#include "../util/problem.hpp"
#include "../util/timer.hpp"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_circle.cpp"
#include "solve_init.cpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  CN_FR,      // Init:CoordinateNewton / Optimize:FR
  CN_L_BFGS,  // Init:CoordinateNewton / Optimize:L_BFGS
  SA_FR,      // Init:SAInitialization / Optimize:FR
  SA_L_BFGS   // Init:SAInitialization / Optimize:L_BFGS
};

std::string MethodStr[6] = {"FR", "L-BFGS", "CN-FR", "CN-L-BFGS", "SA-FR", "SA-L-BFGS"};

std::pair<std::vector<std::pair<double, double>>, std::vector<Eigen::VectorXf>> solve(
    const Method method, const Problem& problem, const bool measureTime, const int seed,
    const int MAX_ITER) {
  Timer timer;
  std::vector<Eigen::VectorXf> positions;

  if (method == CN_FR || method == CN_L_BFGS) {
    solve_init(problem, measureTime, seed, positions, timer);
  } else if (method == SA_FR || method == SA_L_BFGS) {
    solve_circle(problem, seed, positions, timer);
  } else {
    timer.start();
    std::srand(seed);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < int(position.size()); i++) position[i] = std::abs(position[i]);
    positions.push_back(position);
    timer.stop();
  }
  assert(!measureTime || positions.size() == 1);
  std::vector<std::pair<double, double>> hist;
  hist.emplace_back(problem.calcScore(positions.back(), true), timer.sec());

  if (method == L_BFGS || method == CN_L_BFGS || method == SA_L_BFGS) {
    solve_LBFGS<FunctionFR>(problem, positions, hist, timer, MAX_ITER);
  } else {
    solve_FR(problem, positions, hist, timer, MAX_ITER);
  }
  assert(!timer.isMeasuring);
  assert(hist.size() == positions.size());

  return {hist, positions};
}
