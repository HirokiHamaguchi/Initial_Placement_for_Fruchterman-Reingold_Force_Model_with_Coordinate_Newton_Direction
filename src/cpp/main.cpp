#include "include/LBFGS.h"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_RS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

enum Method {
  FR,         // Fruchterman-Reingold method
  L_BFGS,     // Limited-memory Broyden-Fletcher-Goldfarb-Shanno method
  KK_RS,      // Init:Kamada-Kawai / Optimize:Random Subspace Newton
  RS_FR,      // Init:RS / Optimize:FR
  RS_L_BFGS,  // Init:RS / Optimize:L_BFGS
};

int main() {
  // Step1 : Read input

  // const Method method = RS_L_BFGS;
  // const Method method = L_BFGS;
  const Method method = KK_RS;

  // Problem problem("494_bus");
  // Problem problem("dwt_162");
  Problem problem("balanced_tree_2_7");

  bool measureTime = false;

  // Step2 : Solve by each method
  std::vector<Eigen::VectorXf> positions;

  auto t0 = std::chrono::high_resolution_clock::now();

  if (method == RS_FR || method == RS_L_BFGS) {
    positions = solve_init(problem, measureTime);
  } else {
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) position[i] = std::abs(position[i]);
    positions.push_back(position);
  }
  if (method == KK_RS) {
    positions = solve_LBFGS<FunctionKK>(problem, positions.back(), measureTime);
  }

  std::vector<Eigen::VectorXf> positions2;
  if (method == L_BFGS || method == RS_L_BFGS) {
    positions2 = solve_LBFGS<FunctionFR>(problem, positions.back(), measureTime);
  } else if (method == FR || method == RS_FR) {
    positions2 = solve_FR(problem, positions.back(), measureTime);
  } else if (method == KK_RS) {
    positions.push_back(positions.back());
    problem.optimalScaling(positions.back());
    positions2 = solve_RS(problem, positions.back(), measureTime);
    // positions2 = solve_LBFGS<FunctionFR>(problem, positions.back(), measureTime);
  }
  positions.insert(positions.end(), positions2.begin(), positions2.end());

  auto t1 = std::chrono::high_resolution_clock::now();

  // Step3 : Output
  std::cout << "Total Elapsed time: " << std::chrono::duration<double>(t1 - t0).count()
            << " seconds\n"
            << std::endl;

  std::cout << "Final Score: " << problem.calcScore(positions.back()) << std::endl;

  problem.printOutput(positions);

  return 0;
}
