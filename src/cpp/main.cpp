#include "include/LBFGS.h"
#include "solve_FR.cpp"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

enum Method {
  FR,
  L_BFGS,
  RS_FR,
  RS_L_BFGS,
};

int main() {
  // Step1 : Read input
  const Method method = RS_FR;
  // const Method method = FR;
  bool measureTime = false;
  // Problem problem("494_bus");
  // Problem problem("dwt_162");
  Problem problem("balanced_tree_2_5");

  // Step2 : Solve by each method
  std::vector<Eigen::VectorXf> positions;

  auto t0 = std::chrono::high_resolution_clock::now();

  if (method == RS_L_BFGS || method == RS_FR) {
    positions = solve_init(problem, measureTime);
  } else {
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) position[i] = std::abs(position[i]);
    positions.push_back(position);
  }

  if (method == L_BFGS || method == RS_L_BFGS) {
    auto positions2 = solve_LBFGS(problem, positions.back(), measureTime);
    positions.insert(positions.end(), positions2.begin(), positions2.end());
  } else {
    auto positions2 = solve_FR(problem, positions.back(), 100, 1e-3);
    positions.insert(positions.end(), positions2.begin(), positions2.end());
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  // Step3 : Output
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
                   1000.0
            << " seconds\n"
            << std::endl;

  problem.printOutput(positions);

  return 0;
}
