#include "../../../src/cpp/solve/solve.cpp"

int main() {
  const int MAX_ITER = 100;
  std::string matrixName = "jagmesh1";
  Problem problem(matrixName);
  Method method = CN_L_BFGS;
  int seed = 1;

  auto [hist, positions] = solve(method, problem, false, seed, MAX_ITER);
  problem.printOutput(positions, "result.out");

  return 0;
}
