#include "include/LBFGS.h"
#include "solve_LBFGS.cpp"
#include "solve_init.cpp"
#include "util/function.hpp"
#include "util/problem.hpp"

int main() {
  // Step1 : Read input
  Problem problem("jagmesh1");

  // Step2 : Solve
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
    positions = solve_LBFGS(problem, position);
  } else {
    // Step2.1 : Initialization by RandomSubspace
    positions = solve_init(problem);
    // Step2.2 : L-BFGS
    // auto positions2 = solve_LBFGS(problem, positions[positions.size() - 1]);
    // positions.insert(positions.end(), positions2.begin(), positions2.end());
  }

  // Step3 : Output
  problem.printOutput(positions);

  return 0;
}
