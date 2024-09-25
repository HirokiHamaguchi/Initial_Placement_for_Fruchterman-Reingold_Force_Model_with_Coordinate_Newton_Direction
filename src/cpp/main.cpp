#include "include/LBFGS.h"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

int main() {
  // Step1 : Read input
  Problem problem("jagmesh1");

  // Step2 : Solve
  Eigen::VectorXf positions;
  if (false) {
    // Step2.1 : Initialization by Random
    std::srand(0);
    positions = Eigen::VectorXf::Random(2 * problem.n);
    for (int i = 0; i < positions.size(); ++i) {
      assert(-1 <= positions[i] && positions[i] <= 1);
      positions[i] = std::abs(positions[i]);
    }
    // Step2.2 : L-BFGS
    positions = solve_LBFGS(problem, positions);
  } else {
    // Step2.1 : Initialization by RandomSubspace
    positions = solve_init(problem);
    // Step2.2 : L-BFGS
    // positions = solve_LBFGS(problem, positions);
  }

  // Step3 : Output
  problem.printOutput(positions);

  return 0;
}
