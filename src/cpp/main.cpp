#include "solve/solve.cpp"

int main() {
  Method method = CN_L_BFGS;
  Problem problem("jagmesh1");
  bool measureTime = false;

  auto [hist, positions] = solve(method, problem, measureTime, 1, 50);
  assert(!hist.empty());
  std::cout << "Elapsed time: " << hist.back().second << " seconds" << std::endl;
  std::cout << "Score: " << hist.back().first << std::endl;
  problem.printOutput(
      positions, "../../out/" + problem.matrixName + "_" + MethodStr[method] + ".out");

  return 0;
}
