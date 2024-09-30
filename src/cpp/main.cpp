#include "include/LBFGS.h"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

int main() {
  // Step1 : Read input
  bool measureTime = false;
  // Problem problem("494_bus");
  Problem problem("dwt_162");

  // Step2 : Solve
  auto t0 = std::chrono::high_resolution_clock::now();
  std::vector<Eigen::VectorXf> positions;
  if (false) {
    // Step2.1 : Initialization by Random
    std::srand(0);
    Eigen::VectorXf position = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < position.size(); ++i) {
      assert(-1.0 <= position[i] && position[i] <= 1.0);
      position[i] = std::abs(position[i]);
    }
    // Step2.2 : L-BFGS
    positions = solve_LBFGS(problem, position, measureTime);
  } else {
    // Step2.1 : Initialization by RandomSubspace
    positions = solve_init(problem, measureTime);
    // Step2.2 : L-BFGS
    auto positions2 =
        solve_LBFGS(problem, positions[positions.size() - 1], measureTime);
    positions.insert(positions.end(), positions2.begin(), positions2.end());
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() /
                   1000.0
            << " seconds\n"
            << std::endl;

  // Step3 : Output
  problem.printOutput(positions);

  return 0;
}
